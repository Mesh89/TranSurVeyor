#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <unistd.h>
#include <climits>
#include <cassert>

#include <htslib/sam.h>
#include <htslib/kseq.h>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "config.h"
#include "cluster.h"
#include "libs/cptl_stl.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"

config_t config;
std::string workdir;
std::mutex mtx;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;
std::unordered_map<std::string, std::pair<char*, size_t> > chrs;

std::ofstream predictions_writer;

const int SMALL_SAMPLE_SIZE = 15;
const int CLUSTER_CANDIDATES = 3;
const double BASE_ACCEPTANCE_THRESHOLD = 0.95;

const int SKIP_READ = -1;


struct clip_cluster_t {
    cluster_t* c;
    std::string seq;

    clip_cluster_t(cluster_t* c, std::string& seq) : c(c), seq(seq) {}
};

struct reads_cluster_t {
    std::vector<bam1_t*> reads;
    std::vector<clip_cluster_t*> clip_clusters;
};


struct region_t {
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    int score = 0;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}
};

struct cc_v_distance_t {
    reads_cluster_t c1, c2;
    int distance;

    cc_v_distance_t(reads_cluster_t c1, reads_cluster_t c2, int distance) : c1(c1), c2(c2), distance(distance) {}
};
bool operator < (const cc_v_distance_t& ccd1, const cc_v_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance < ccd2.distance;
}

void remap_supp(int contig_id, std::unordered_map<std::string, std::string> &mateseqs, std::string &qname,
                const std::string &l_dc_fname);

std::string print_region(region_t region) {
    std::stringstream ss;
    ss << contig_id2name[region.contig_id] << ":" << region.start << "-" << region.end;
    return ss.str();
}

region_t get_region(std::vector<bam1_t*> subcluster, std::string& m_contig_name) {
    bam1_t* leftmost_reverse_mate = NULL, * rightmost_forward_mate = NULL;
    /* gets two regions
     * 1. from left-most reverse read, max_is to the left and max_tra_size to the right
     * 2. from right-most forward read, max_is to the right and max_tra_size to the left
     */
    for (bam1_t* r : subcluster) { // get leftmost reverse read
        if (bam_is_mrev(r) && (leftmost_reverse_mate == NULL || leftmost_reverse_mate->core.mpos > r->core.mpos)) {
            leftmost_reverse_mate = r;
        }
    }
    for (bam1_t* r : subcluster) { // get rightmost forward read
        if (!bam_is_mrev(r) && (rightmost_forward_mate == NULL || rightmost_forward_mate->core.mpos < r->core.mpos)) {
            rightmost_forward_mate = r;
        }
    }

    int start = INT_MAX;
    int end = 0;
    if (leftmost_reverse_mate != NULL) {
        start = std::min(start, leftmost_reverse_mate->core.mpos-config.max_is);
        end = std::max(end, leftmost_reverse_mate->core.mpos+config.max_tra_size);
    }
    if (rightmost_forward_mate != NULL) {
        start = std::min(start, rightmost_forward_mate->core.mpos-config.max_tra_size);
        end = std::max(end, rightmost_forward_mate->core.mpos+config.max_is);
    }

    std::pair<char *, size_t> chr = chrs[m_contig_name];
    int contig_len = chr.second;
    return region_t(contig_name2id[m_contig_name], subcluster[0]->core.mtid, std::max(0,start), std::min(end,contig_len));
}

char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
};

int compute_score_supp(region_t& region, std::vector<bam1_t*>& reads, std::unordered_map<std::string, std::string>& mateseqs,
                   std::vector<clip_cluster_t*> clip_clusters, std::vector<int>* offsets, std::vector<std::string>* cigars,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter, bool do_rc) {

    std::vector<std::string> mates;
    for (bam1_t* r : reads) {
        std::string qname = bam_get_qname(r);
        if (is_samechr(r)) {
            if (r->core.isize > 0) qname += "_2";
            else qname += "_1";
        }
        std::string s = mateseqs[qname];

        // mateseqs contains seqs converted to positive strand
        //  do_rc: pos stable strand = pos unstable strand and neg stable strand = neg unstable strand
        // !do_rc: pos stable strand = neg unstable strand and neg stable strand = pos unstable strand
        if ((do_rc && bam_is_rev(r)) ||      //  do_rc: rc mate when neg stable
            (!do_rc && !bam_is_rev(r))) {    // !do_rc: rc mate when pos strand
            rc(s);
        }

        mates.push_back(s);
    }
    for (clip_cluster_t* cc : clip_clusters) {
        std::string s = cc->seq;
        if (do_rc) {
            rc(s);
        }
        mates.push_back(s);
    }

    int score = 0;
    for (std::string& s : mates) {
        StripedSmithWaterman::Alignment alignment;
        int mask_len = s.length()/2;
        if (mask_len < 15) mask_len = 15;

        aligner.Align(s.c_str(), chrs[contig_id2name[region.contig_id]].first+region.start, region.end-region.start,
                      filter, &alignment, mask_len);
        bool accepted = alignment.sw_score >= s.length();
        accepted &= !is_poly_ACGT(s.c_str()+alignment.query_begin, alignment.query_end-alignment.query_begin+1);
        if (accepted) {
            score += alignment.sw_score;
        }

        if (offsets != NULL) {
            if (accepted) {
                offsets->push_back(alignment.ref_begin);
            } else {
                offsets->push_back(SKIP_READ);
            }
        }
        if (cigars != NULL) {
            // wrapper that returns M in case of = or X
            std::stringstream ss;
            char op = ' '; int len = 0;
            for (uint32_t c : alignment.cigar) {
                if (op != _cigar_int_to_op(c)) {
                    if (op != ' ') ss << len << op;
                    op = _cigar_int_to_op(c);
                    len = cigar_int_to_len(c);
                } else {
                    len += cigar_int_to_len(c);
                }
            }
            ss << len << op;
            cigars->push_back(ss.str());
        }
    }
    return score;
}

void compute_score(region_t& region, std::vector<bam1_t*>& reads, std::unordered_map<std::string, std::string>& mateseqs,
                   std::vector<clip_cluster_t*> clip_clusters, std::vector<int>* offsets, std::vector<std::string>* cigars,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter, bool& is_rc) {

    int score = compute_score_supp(region, reads, mateseqs, clip_clusters, NULL, NULL, aligner, filter, false);
    int rc_score = compute_score_supp(region, reads, mateseqs, clip_clusters, NULL, NULL, aligner, filter, true);
    region.score = std::max(score, rc_score);
    if (score >= rc_score) {
        is_rc = false;
        if (offsets != NULL) {
            compute_score_supp(region, reads, mateseqs, clip_clusters, offsets, cigars, aligner, filter, false);
        }
    } else {
        is_rc = true;
        if (offsets != NULL) {
            compute_score_supp(region, reads, mateseqs, clip_clusters, offsets, cigars, aligner, filter, true);
        }
    }
}

samFile* open_writer(std::string name, bam_hdr_t* header) {
    samFile* remapped_file = sam_open(name.c_str(), "wb");
    if (remapped_file == NULL) {
        throw "Unable to open " + name;
    }
    if (sam_hdr_write(remapped_file, header) != 0) {
        throw "Could not write file " + std::string(remapped_file->fn);
    }
    return remapped_file;
}

prediction_t make_prediction(reads_cluster_t reads, int contig_id, int mcontig_id) {
    auto clip_score = [] (const std::vector<clip_cluster_t*> clips) {
        int sum = 0;
        for (clip_cluster_t* cc : clips) {
            sum += cc->c->a1.sc_reads;
        }
        return sum;
    };

    auto start_comp = [] (const bam1_t* r1, const bam1_t* r2) { return r1->core.pos < r2->core.pos; };
    auto end_comp = [] (const bam1_t* r1, const bam1_t* r2) { return bam_endpos(r1) < bam_endpos(r2); };
    int start_pos = (*std::min_element(reads.reads.begin(), reads.reads.end(), start_comp))->core.pos;
    int end_pos = bam_endpos(*std::max_element(reads.reads.begin(), reads.reads.end(), end_comp));
    anchor_t a1(bam_is_rev(reads.reads[0]) ? 'L' : 'R', contig_id, start_pos, end_pos, clip_score(reads.clip_clusters));

    auto mstart_comp = [] (const bam1_t* r1, const bam1_t* r2) { return r1->core.mpos < r2->core.mpos; };
    auto mend_comp = [] (const bam1_t* r1, const bam1_t* r2) { return get_mate_endpos(r1) < bam_endpos(r2); };
    int mstart_pos = (*std::min_element(reads.reads.begin(), reads.reads.end(), mstart_comp))->core.mpos;
    int mend_pos = get_mate_endpos(*std::max_element(reads.reads.begin(), reads.reads.end(), mend_comp));
    anchor_t a2(bam_is_mrev(reads.reads[0]) ? 'L' : 'R', mcontig_id, mstart_pos, mend_pos, clip_score(reads.clip_clusters));

    cluster_t* c = new cluster_t(a1, a2, DISC_TYPES.DC, reads.reads.size());
    prediction_t prediction(c, DISC_TYPES.DC);
    delete c;
    return prediction;
}

std::atomic<int> loc_pred_id;

void remap_cluster(reads_cluster_t r_cluster, reads_cluster_t l_cluster, std::vector<bam1_t*>& kept,
                   int contig_id, bam_hdr_t* header, std::unordered_map<std::string, std::string>& mateseqs,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& aligner_to_base) {
    std::vector<region_t> regions;

    // sort by mate chr and mate pos
    std::vector<bam1_t*> full_cluster;
    full_cluster.insert(full_cluster.end(), l_cluster.reads.begin(), l_cluster.reads.end());
    full_cluster.insert(full_cluster.end(), r_cluster.reads.begin(), r_cluster.reads.end());
    sort(full_cluster.begin(), full_cluster.end(), [] (bam1_t* r1, bam1_t* r2) {
        if (r1->core.mtid != r2->core.mtid) return r1->core.mtid < r2->core.mtid;
        else return r1->core.mpos < r2->core.mpos;
    });

    // cluster the reads according to the mates
    std::vector<bam1_t*> subcluster;
    for (bam1_t* r : full_cluster) {
        if (!subcluster.empty() && (subcluster[0]->core.mtid != r->core.mtid ||
                r->core.mpos-subcluster[0]->core.mpos > config.max_is)) {
            std::string m_contig_name = std::string(header->target_name[subcluster[0]->core.mtid]);
            regions.push_back(get_region(subcluster, m_contig_name));
            subcluster.clear();
        }
        subcluster.push_back(r);
    }
    if (!subcluster.empty()) {
        std::string m_contig_name = std::string(header->target_name[subcluster[0]->core.mtid]);
        regions.push_back(get_region(subcluster, m_contig_name));
    }

    StripedSmithWaterman::Filter filter, filter_w_cigar;
    filter_w_cigar.report_cigar = true;
    bool is_rc;

    std::vector<clip_cluster_t*> clip_clusters;
    clip_clusters.insert(clip_clusters.end(), l_cluster.clip_clusters.begin(), l_cluster.clip_clusters.end());
    clip_clusters.insert(clip_clusters.end(), r_cluster.clip_clusters.begin(), r_cluster.clip_clusters.end());

    // if too much regions and too many reads, subsample
    if (full_cluster.size() > SMALL_SAMPLE_SIZE && regions.size() > CLUSTER_CANDIDATES) {
        std::vector<bam1_t*> small_sample(full_cluster);
        std::random_shuffle(small_sample.begin(), small_sample.end());
        if (small_sample.size() > SMALL_SAMPLE_SIZE) {
            small_sample.erase(small_sample.begin() + SMALL_SAMPLE_SIZE, small_sample.end());
        }

        // compute best score
        for (int i = 0; i < regions.size(); i++) {
            compute_score(regions[i], small_sample, mateseqs, clip_clusters, NULL, NULL, aligner, filter, is_rc);
        }
        sort(regions.begin(), regions.end(), [] (region_t r1, region_t r2) {return r1.score > r2.score;});

        regions.erase(regions.begin()+CLUSTER_CANDIDATES, regions.end());
    }

    // compute best score
    for (int i = 0; i < regions.size(); i++) {
        compute_score(regions[i], full_cluster, mateseqs, clip_clusters, NULL, NULL, aligner, filter, is_rc);
    }
    sort(regions.begin(), regions.end(), [] (region_t r1, region_t r2) {return r1.score > r2.score;});

    region_t best_region = regions[0];

    // get base region
    sort(full_cluster.begin(), full_cluster.end(), [] (bam1_t* r1, bam1_t* r2) {
        return r1->core.pos < r2->core.pos;
    });
    int start = full_cluster[0]->core.pos - config.max_is;
    int end = bam_endpos(full_cluster[full_cluster.size()-1]) + config.max_is;
    int contig_len = chrs[contig_id2name[contig_id]].second;
    region_t base_region(contig_id, full_cluster[0]->core.tid, std::max(0,start), std::min(end,contig_len));

    compute_score(base_region, full_cluster, mateseqs, clip_clusters, NULL, NULL, aligner_to_base, filter, is_rc);

    if (base_region.score >= best_region.score*BASE_ACCEPTANCE_THRESHOLD) {
        return;
    }

    reads_cluster_t pos_cluster, neg_cluster;
    std::vector<int> offsets;
    std::vector<std::string> cigars;
    compute_score(best_region, full_cluster, mateseqs, clip_clusters, &offsets, &cigars, aligner, filter_w_cigar, is_rc);
    for (int i = 0; i < full_cluster.size(); i++) {
        if (offsets[i] == SKIP_READ) continue; // TODO: mem leak

        bam1_t* r = full_cluster[i];
        r->core.mtid = best_region.original_bam_id;
        r->core.mpos = best_region.start + offsets[i];
        if (is_rc == bam_is_rev(r)) {
            r->core.flag |= BAM_FMREVERSE; //sets flag to true
            assert(bam_is_mrev(r));
        } else {
            r->core.flag &= ~BAM_FMREVERSE; //sets flag to false
            assert(!bam_is_mrev(r));
        }

        if (!bam_is_rev(r)) {
            pos_cluster.reads.push_back(r);
        } else {
            neg_cluster.reads.push_back(r);
        }

        bam_aux_update_str(r, "MC", cigars[i].length()+1, cigars[i].c_str());

        kept.push_back(r);
    }
    for (int i = 0; i < clip_clusters.size(); i++) {
        if (offsets[full_cluster.size()+i] == SKIP_READ) continue;

        if (clip_clusters[i]->c->a1.dir == 'R') {
            pos_cluster.clip_clusters.push_back(clip_clusters[i]);
        } else {
            neg_cluster.clip_clusters.push_back(clip_clusters[i]);
        }
    }

    int pred_id = loc_pred_id++;
    if (!pos_cluster.reads.empty()) {
        prediction_t pred = make_prediction(pos_cluster, contig_id, best_region.contig_id);
        pred.id = pred_id;
        mtx.lock();
        predictions_writer << pred.to_str() << "\n";
        mtx.unlock();
    }
    if (!neg_cluster.reads.empty()) {
        prediction_t pred = make_prediction(neg_cluster, contig_id, best_region.contig_id);
        pred.id = pred_id;
        mtx.lock();
        predictions_writer << pred.to_str() << "\n";
        mtx.unlock();
    }

    sort(kept.begin(), kept.end(), [] (bam1_t* r1, bam1_t* r2) {return get_endpoint(r1) < get_endpoint(r2);});
}

int find(int* parents, int i) {
    int root = i;
    while (root != parents[root]) {
        root = parents[root];
    }
    while (i != root) {
        int newp = parents[i];
        parents[i] = root;
        i = newp;
    }
    return root;
}
void merge(int* parents, int* sizes, int x, int y) {
    int i = find(parents, x);
    int j = find(parents, y);
    if (i == j) return;

    if (sizes[i] < sizes[j]) {
        parents[i] = j;
        sizes[j] += sizes[i];
    } else {
        parents[j] = i;
        sizes[i] += sizes[j];
    }
}

void remove_cluster_from_mm(std::multimap<int, cluster_t*>& mm, cluster_t* c, int pos) {
    auto bounds = mm.equal_range(pos);
    for (auto it = bounds.first; it != bounds.second; it++) {
        if (it->second == c) {
            mm.erase(it);
            break;
        }
    }
}
void remove_cluster_from_mm(std::multimap<int, cluster_t*>& mm, cluster_t* c) {
    remove_cluster_from_mm(mm, c, c->a1.start);
    remove_cluster_from_mm(mm, c, c->a1.end);
}

std::vector<reads_cluster_t> cluster_reads(open_samFile_t* dc_file, int contig_id,
                                           std::unordered_map<std::string, std::string>& mateseqs,
                                           std::vector<clip_cluster_t*>& clip_clusters) {

    std::string contig = contig_id2name[contig_id];
    hts_itr_t* iter = sam_itr_querys(dc_file->idx, dc_file->header, contig.c_str());
    bam1_t* read = bam_init1();

    std::vector<cluster_t*> clusters;
    std::multimap<int, cluster_t*> clusters_map;
    std::vector<bam1_t*> reads;
    while (sam_itr_next(dc_file->file, iter, read) >= 0) {
        std::string qname = bam_get_qname(read);
        if (mateseqs.count(qname) == 0 && mateseqs.count(qname + "_1") == 0 &&
            mateseqs.count(qname + "_2") == 0) continue; // mateseq not present

        std::string mate_read = mateseqs[qname];
        if (is_poly_ACGT(read) || is_poly_ACGT(mate_read.c_str(), mate_read.length())) continue;

        bool rev = bam_is_rev(read);
        anchor_t a(rev ? 'L' : 'R', contig_id, read->core.pos, bam_endpos(read), 0);
        cluster_t* c = new cluster_t(a, a, DISC_TYPES.DC, 1);

        c->id = clusters.size();
        clusters.push_back(c);
        reads.push_back(bam_dup1(read));
    }
    for (int i = 0; i < clip_clusters.size(); i++) {
        cluster_t* c = clip_clusters[i]->c;
        c->id = i + reads.size(); // clips have negative ids
        clusters.push_back(c);
    }

    sam_itr_destroy(iter);
    bam_destroy1(read);

    if (clusters.empty()) return std::vector<reads_cluster_t>();

    // union-find datastructure
    int n_reads = reads.size() + clip_clusters.size();
    int* parents = new int[n_reads], * sizes = new int[n_reads];
    for (int i = 0; i < n_reads; i++) {
        parents[i] = i;
        sizes[i] = 1;
    }

    sort(clusters.begin(), clusters.end(), [](const cluster_t* c1, const cluster_t* c2) {
        return c1->a1 < c2->a1;
    });

    // merge first equal clusters
    int prev;
    do {
        prev = clusters.size();
        for (int i = 0; i < clusters.size()-1; i++) {
            cluster_t* c1 = clusters[i], * c2 = clusters[i+1];
            if (c1 != NULL && c1->a1.start == c2->a1.start && c1->a1.end == c2->a1.end) {
                cluster_t* new_cluster = cluster_t::merge(c1, c2);
                new_cluster->id = std::min(c1->id, c2->id);
                assert(std::min(c1->id, c2->id) >= 0);
                merge(parents, sizes, c1->id, c2->id);
                clusters[i] = new_cluster;
                clusters[i+1] = NULL;
            }
        }

        clusters.erase(std::remove(clusters.begin(), clusters.end(), (cluster_t*) NULL), clusters.end());
    } while (prev != clusters.size());

    for (cluster_t* c : clusters) {
        clusters_map.insert(std::make_pair(c->a1.start, c));
        clusters_map.insert(std::make_pair(c->a1.end, c));
    }


    std::vector<int> max_dists;
//    for (int i = -95; i < 100; i+=5) max_dists.push_back(i);
//    for (int i = 100; i < config.max_is; i+=100) max_dists.push_back(i);
    max_dists.push_back(config.max_is);

    for (int max_dist : max_dists) {
        std::priority_queue<cc_distance_t> pq;
        for (cluster_t* c1 : clusters) {
            if (c1->dead) continue;
            auto end = clusters_map.upper_bound(c1->a1.end+std::max(max_dist, 0));
            for (auto map_it = clusters_map.lower_bound(c1->a1.start); map_it != end; map_it++) {
                cluster_t* c2 = map_it->second;
                int dist = anchor_t::distance(c1->a1, c2->a1);
                if (c1 != c2 && anchor_t::can_merge(c1->a1, c2->a1, config) && c1->a1.start <= c2->a1.start && dist <= max_dist) {
                    pq.push(cc_distance_t(dist, c1, c2));
                }
            }
        }

        while (!pq.empty()) {
            cc_distance_t ccd = pq.top();
            pq.pop();

            if (ccd.c1->dead || ccd.c2->dead) continue;

            cluster_t* new_cluster = cluster_t::merge(ccd.c1, ccd.c2);
            new_cluster->id = std::min(ccd.c1->id, ccd.c2->id); // clip clusters have negative id
            assert(std::min(ccd.c1->id, ccd.c2->id) >= 0);
            merge(parents, sizes, ccd.c1->id, ccd.c2->id);
            clusters.push_back(new_cluster);

            ccd.c1->dead = true;
            remove_cluster_from_mm(clusters_map, ccd.c1);
            ccd.c2->dead = true;
            remove_cluster_from_mm(clusters_map, ccd.c2);

            auto end = clusters_map.upper_bound(new_cluster->a1.end + max_dist);
            for (auto map_it = clusters_map.lower_bound(new_cluster->a1.start - max_dist);
                 map_it != end; map_it++) {
                if (!map_it->second->dead && cluster_t::can_merge(new_cluster, map_it->second, config)) {
                    pq.push(cc_distance_t(cluster_t::distance(new_cluster, map_it->second), new_cluster,
                                          map_it->second));
                }
            }
            clusters_map.insert(std::make_pair(new_cluster->a1.start, new_cluster));
            clusters_map.insert(std::make_pair(new_cluster->a1.end, new_cluster));
        }
    }

    // for each set of reads, make a cluster-vector
    std::vector<reads_cluster_t> read_clusters;
    for (int i = 0; i < n_reads; i++) {
        read_clusters.push_back(reads_cluster_t());
    }
    for (int i = 0; i < reads.size(); i++) {
        read_clusters[find(parents, i)].reads.push_back(reads[i]);
    }
    for (int i = 0; i < clip_clusters.size(); i++) {
        read_clusters[find(parents, i+reads.size())].clip_clusters.push_back(clip_clusters[i]);
    }

    // remove clusters of size < 2 (with no mem leaks)
    for (reads_cluster_t& rc : read_clusters) {
        if (rc.reads.size() == 1) bam_destroy1(rc.reads[0]);
    }
    read_clusters.erase(std::remove_if(read_clusters.begin(), read_clusters.end(),
                                       [](reads_cluster_t rc) {return rc.reads.size() <= 1;}), read_clusters.end());

    for (cluster_t* c : clusters) delete c;

    delete[] parents;
    delete[] sizes;

    return read_clusters;
}

void write_and_index_file(std::vector<bam1_t*>& reads, std::string fname, bam_hdr_t* header) {
    samFile* file = open_writer(fname, header);
    if (file == NULL) {
        throw "Unable to open " + fname;
    }

    // write reads
    sort(reads.begin(), reads.end(), [](bam1_t *r1, bam1_t *r2) { return r1->core.pos < r2->core.pos; });
    for (bam1_t* r : reads) {
        int ok = sam_write1(file, header, r);
        if (ok < 0) throw "Unable to write to " + fname;
    }

    sam_close(file);

    file = sam_open(fname.c_str(), "r");

    int code = sam_index_build(fname.c_str(), 0);
    if (code != 0) {
        throw "Cannot index " + fname;
    }

    sam_close(file);
}

void remap(int id, int contig_id,   std::vector<clip_cluster_t*>& r_clip_clusters,
                                    std::vector<clip_cluster_t*>& l_clip_clusters) {
    mtx.lock();
    std::cout << "Remapping DC for " << contig_id << " (" << contig_id2name[contig_id] << ")" << std::endl;
    mtx.unlock();

    StripedSmithWaterman::Aligner aligner(2, 2, 3, 1, false);
    StripedSmithWaterman::Aligner aligner_to_base(2, 2, 3, 1, true);

    std::unordered_map<std::string, std::string> mateseqs;
    std::ifstream mateseqs_fin(workdir + "/workspace/" + std::to_string(contig_id) + "-MATESEQS");
    std::string qname, seq;
    while (mateseqs_fin >> qname >> seq) {
        mateseqs[qname] = seq;
    }
    mateseqs_fin.close();

    std::string l_dc_fname = workdir + "/workspace/L" + std::to_string(contig_id) + "-DC.noremap.bam";
    std::string r_dc_fname = workdir + "/workspace/R" + std::to_string(contig_id) + "-DC.noremap.bam";
    open_samFile_t* l_dc_file = open_samFile(l_dc_fname.c_str(), true);
    open_samFile_t* r_dc_file = open_samFile(r_dc_fname.c_str(), true);
    std::vector<reads_cluster_t> l_clusters = cluster_reads(l_dc_file, contig_id, mateseqs, l_clip_clusters);
    std::vector<reads_cluster_t> r_clusters = cluster_reads(r_dc_file, contig_id, mateseqs, r_clip_clusters);

    std::vector<bam1_t*> l_reads_to_write, r_reads_to_write;

    std::priority_queue<cc_v_distance_t> pq;
    auto score_f = [](const reads_cluster_t v1, const reads_cluster_t v2) {return v1.reads.size()*v2.reads.size();};

    std::set<int> r_clusters_available, l_clusters_available;
    std::multimap<int, reads_cluster_t> l_clusters_map;
    for (reads_cluster_t& l_cluster : l_clusters) {
        int pos = l_cluster.reads[0]->core.pos;
        l_clusters_map.insert(std::make_pair(pos, l_cluster));
        l_clusters_available.insert(pos);
    }
    for (reads_cluster_t& r_cluster : r_clusters) {
        int pos = bam_endpos(r_cluster.reads[r_cluster.reads.size()-1]);
        auto begin = l_clusters_map.lower_bound(pos-50);
        auto end = l_clusters_map.upper_bound(pos+config.max_is);

        reads_cluster_t* l_cluster = NULL;
        for (auto it = begin; it != end; it++) {
            if (l_cluster == NULL || it->second.reads.size() > l_cluster->reads.size()) {
                l_cluster = &it->second;
            }
        }
        if (l_cluster == NULL) continue;

        r_clusters_available.insert(pos);

        pq.push(cc_v_distance_t(r_cluster, *l_cluster, score_f(r_cluster, *l_cluster)));
    }

    while (!pq.empty()) {
        cc_v_distance_t cc_v_distance = pq.top();
        pq.pop();

        reads_cluster_t c1 = cc_v_distance.c1;
        reads_cluster_t c2 = cc_v_distance.c2;
        int r_pos = bam_endpos(*(c1.reads.rbegin()));
        int l_pos = c2.reads[0]->core.pos;

        auto r_it = r_clusters_available.find(r_pos);
        auto l_it = l_clusters_available.find(l_pos);
        if (r_it == r_clusters_available.end() || l_it == l_clusters_available.end()) continue;

        r_clusters_available.erase(r_it);
        l_clusters_available.erase(l_it);

        // remap clusters
        std::vector<bam1_t*> to_write;
        remap_cluster(c1, c2, to_write, contig_id, r_dc_file->header, mateseqs, aligner, aligner_to_base);

        for (bam1_t* r : to_write) {
            if (bam_is_rev(r)) {
                l_reads_to_write.push_back(r);
            } else {
                r_reads_to_write.push_back(r);
            }
        }
    }


    std::string l_dc_remapped_fname = workdir + "/workspace/L" + std::to_string(contig_id) + "-DC.remap.bam";
    write_and_index_file(l_reads_to_write, l_dc_remapped_fname, l_dc_file->header);

    std::string r_dc_remapped_fname = workdir + "/workspace/R" + std::to_string(contig_id) + "-DC.remap.bam";
    write_and_index_file(r_reads_to_write, r_dc_remapped_fname, r_dc_file->header);

    close_samFile(l_dc_file);
    close_samFile(r_dc_file);

    // destroy reads
    for (reads_cluster_t l_cluster : l_clusters) {
        for (bam1_t* r : l_cluster.reads) {
            bam_destroy1(r);
        }
    }
    for (reads_cluster_t r_cluster : r_clusters) {
        for (bam1_t* r : r_cluster.reads) {
            bam_destroy1(r);
        }
    }
}

int main(int argc, char* argv[]) {

    workdir = std::string(argv[1]);
    std::string workspace = workdir + "/workspace";
    std::string reference_fname  = std::string(argv[2]);

    FILE* fastaf = fopen(reference_fname.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(fastaf));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        chrs[std::string(seq->name.s)] = std::make_pair(new char[seq->seq.l+1], seq->seq.l);
        strcpy(chrs[std::string(seq->name.s)].first, seq->seq.s);
    }
    kseq_destroy(seq);

    config = parse_config(workdir + "/config.txt");

    // we explicitly store the contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }

    predictions_writer.open(workspace + "/predictions.raw");

    std::string clip_name, clip_seq;
    std::ifstream clips_fin(workspace + "/CLIPS.fa");
    std::unordered_map<int, std::vector<clip_cluster_t*> > r_clip_clusters, l_clip_clusters;
    while (getline(clips_fin, clip_name)) {
        getline(clips_fin, clip_seq);
        int anchor_contig_id, anchor_start, anchor_end; char anchor_dir; int anchor_sc_reads;
        sscanf(clip_name.c_str()+1, "%d_%d_%d_%c_%d", &anchor_contig_id, &anchor_start, &anchor_end, &anchor_dir, &anchor_sc_reads);
        anchor_t a(anchor_dir, anchor_contig_id, anchor_start, anchor_end, anchor_sc_reads);
        cluster_t* c = new cluster_t(a, a, DISC_TYPES.DC, 0);
        if (anchor_dir == 'L') {
            l_clip_clusters[anchor_contig_id].push_back(new clip_cluster_t(c, clip_seq));
        } else {
            r_clip_clusters[anchor_contig_id].push_back(new clip_cluster_t(c, clip_seq));
        }
    }

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (contig_id = 1; contig_id < contig_id2name.size(); contig_id++) {
        std::future<void> future = thread_pool.push(remap, contig_id, r_clip_clusters[contig_id], l_clip_clusters[contig_id]);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    // free memory
    for (auto& e : l_clip_clusters) {
        for (clip_cluster_t* cc : e.second) {
            delete cc;
        }
    }
    for (auto& e : r_clip_clusters) {
        for (clip_cluster_t* cc : e.second) {
            delete cc;
        }
    }

    for (auto& chr : chrs) {
        delete[] chr.second.first;
    }

    predictions_writer.close();
}