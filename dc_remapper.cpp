#include <iostream>
#include <unordered_map>
#include <map>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <cassert>
#include <htslib/hts.h>

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

const int SMALL_SAMPLE_SIZE = 15;
const int CLUSTER_CANDIDATES = 3;
const double BASE_ACCEPTANCE_THRESHOLD = 0.85;

const int SKIP_READ = -1;

struct region_t {
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    int score = 0;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}
};

std::string print_region(region_t region) {
    std::stringstream ss;
    ss << contig_id2name[region.contig_id] << ":" << region.start << "-" << region.end;
    return ss.str();
}

region_t get_region(std::vector<bam1_t*> subcluster, std::string& m_contig_name) {
    int start = subcluster[0]->core.mpos - config.max_is;
    int end = get_mate_endpos(subcluster[subcluster.size()-1]) + config.max_is;
    std::pair<char *, size_t> chr = chrs[m_contig_name];
    int contig_len = chr.second;
    return region_t(contig_name2id[m_contig_name], subcluster[0]->core.mtid, std::max(0,start), std::min(end,contig_len));
}

char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
};

int compute_score_supp(region_t& region, std::vector<bam1_t*>& reads, std::unordered_map<std::string, std::string>& mateseqs,
                   std::vector<int>* offsets, std::vector<std::string>* cigars,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter, bool do_rc) {
    int score = 0;
    for (bam1_t* r : reads) {
        std::string qname = bam_get_qname(r);
        if (is_samechr(r)) {
            if (r->core.isize > 0) qname += "_2";
            else qname += "_1";
        }
        std::string s = mateseqs[qname];
        StripedSmithWaterman::Alignment alignment;
        int mask_len = s.length()/2;
        if (mask_len < 15) mask_len = 15;

        if (do_rc) {
            rc(s);
        }
        aligner.Align(s.c_str(), chrs[contig_id2name[region.contig_id]].first+region.start, region.end-region.start,
                      filter, &alignment, mask_len);
        bool accepted = alignment.sw_score >= r->core.l_qseq;
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
                   std::vector<int>* offsets, std::vector<std::string>* cigars,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter, bool& is_rc) {
    int score = compute_score_supp(region, reads, mateseqs, NULL, NULL, aligner, filter, false);
    int rc_score = compute_score_supp(region, reads, mateseqs, NULL, NULL, aligner, filter, true);
    region.score = std::max(score, rc_score);
    if (score >= rc_score) {
        is_rc = false;
        if (offsets != NULL) {
            compute_score_supp(region, reads, mateseqs, offsets, cigars, aligner, filter, false);
        }
    } else {
        is_rc = true;
        if (offsets != NULL) {
            compute_score_supp(region, reads, mateseqs, offsets, cigars, aligner, filter, true);
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

void remap_cluster(std::vector<bam1_t*>& cluster, int contig_id, bam_hdr_t* header, bool rev,
                   std::unordered_map<std::string, std::string>& mateseqs,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& aligner_to_base,
                   samFile* dc_remapped_file) {
    if (cluster.size() == 1) return;

    sort(cluster.begin(), cluster.end(), [] (bam1_t* r1, bam1_t* r2) {
        if (r1->core.mtid != r2->core.mtid) return r1->core.mtid < r2->core.mtid;
        else return r1->core.mpos < r2->core.mpos;
    });

    std::vector<region_t> regions;

    std::vector<bam1_t*> subcluster;
    for (bam1_t* r : cluster) {
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

    if (regions.empty()) return;

    StripedSmithWaterman::Filter filter, filter_w_cigar;
    filter_w_cigar.report_cigar = true;
    bool is_rc;

    // TODO: make N match
    if (cluster.size() > SMALL_SAMPLE_SIZE && regions.size() > CLUSTER_CANDIDATES) {
        std::vector<bam1_t*> small_sample(cluster);
        std::random_shuffle(small_sample.begin(), small_sample.end());
        if (small_sample.size() > SMALL_SAMPLE_SIZE) {
            small_sample.erase(small_sample.begin() + SMALL_SAMPLE_SIZE, small_sample.end());
        }

        // compute best score
        for (int i = 0; i < regions.size(); i++) {
            compute_score(regions[i], small_sample, mateseqs, NULL, NULL, aligner, filter, is_rc);
        }
        sort(regions.begin(), regions.end(), [] (region_t r1, region_t r2) {return r1.score > r2.score;});

        regions.erase(regions.begin()+CLUSTER_CANDIDATES, regions.end());
    }

    // compute best score
    for (int i = 0; i < regions.size(); i++) {
        compute_score(regions[i], cluster, mateseqs, NULL, NULL, aligner, filter, is_rc);
    }
    sort(regions.begin(), regions.end(), [] (region_t r1, region_t r2) {return r1.score > r2.score;});

    region_t best_region = regions[0];

    // get base region
    int start = cluster[0]->core.pos - config.max_is;
    int end = bam_endpos(cluster[cluster.size()-1]) + config.max_is;
    int contig_len = chrs[contig_id2name[contig_id]].second;
    region_t base_region(contig_id, cluster[0]->core.tid, std::max(0,start), std::min(end,contig_len));

    compute_score(base_region, cluster, mateseqs, NULL, NULL, aligner_to_base, filter, is_rc);

    if (base_region.score >= best_region.score*BASE_ACCEPTANCE_THRESHOLD) {
        cluster.clear();
        return;
    }

    std::vector<bam1_t*> kept;
    std::vector<int> offsets;
    std::vector<std::string> cigars;
    compute_score(best_region, cluster, mateseqs, &offsets, &cigars, aligner, filter_w_cigar, is_rc);
    for (int i = 0; i < cluster.size(); i++) {
        if (offsets[i] == SKIP_READ) continue;

        bam1_t* r = cluster[i];
        r->core.mtid = best_region.original_bam_id;
        r->core.mpos = best_region.start + offsets[i];
        if (is_rc) {
            r->core.flag |= BAM_FMREVERSE; //sets flag to true
                assert(bam_is_mrev(r));
        } else {
            r->core.flag &= ~BAM_FMREVERSE; //sets flag to false
            assert(!bam_is_mrev(r));
        }

        std::string dir = rev ? "L" : "R";
        bam_aux_update_str(r, "MC", cigars[i].length()+1, cigars[i].c_str());

        kept.push_back(r);
    }

    cluster.swap(kept);
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

void remap(int id, int contig_id, bool rev) {
    mtx.lock();
    std::cout << "Remapping DC for " << contig_id << " (" << contig_id2name[contig_id] << ")" << std::endl;
    mtx.unlock();

    std::unordered_map<std::string, std::string> mateseqs;
    std::ifstream mateseqs_fin(workdir + "/workspace/" + std::to_string(contig_id) + "-MATESEQS");
    std::string qname, seq;
    while (mateseqs_fin >> qname >> seq) {
        mateseqs[qname] = seq;
    }
    mateseqs_fin.close();

    std::string dir = rev ? "L" : "R";
    std::string dc_fname = workdir + "/workspace/" + dir + std::to_string(contig_id) + "-DC.noremap.bam";
    samFile* dc_file = sam_open(dc_fname.c_str(), "r");
    if (dc_file == NULL) {
        throw "Unable to open " + dc_fname;
    }

    int code = sam_index_build(dc_fname.c_str(), 0);
    if (code != 0) {
        throw "Cannot index " + dc_fname;
    }

    hts_idx_t* idx = sam_index_load(dc_file, dc_fname.c_str());
    if (idx == NULL) {
        throw "Unable to open index for " + dc_fname;
    }

    bam_hdr_t* header = sam_hdr_read(dc_file);
    if (header == NULL) {
        throw "Unable to open header for " + dc_fname;
    }

    std::string dc_remapped_fname = workdir + "/workspace/" + dir + std::to_string(contig_id) + "-DC.remap.bam";
    samFile* dc_remapped_file = open_writer(dc_remapped_fname, header);
    if (dc_remapped_file == NULL) {
        throw "Unable to open " + dc_remapped_fname;
    }

    StripedSmithWaterman::Aligner aligner(2,2,3,1,false);
    StripedSmithWaterman::Aligner aligner_to_base(2,2,3,1,true);

    std::string contig = contig_id2name[contig_id];
    hts_itr_t* iter = sam_itr_querys(idx, header, contig.c_str());
    bam1_t* read = bam_init1();

    std::vector<cluster_t*> clusters;
    std::multimap<int, cluster_t*> clusters_map;
    std::vector<bam1_t*> reads;
    while (sam_itr_next(dc_file, iter, read) >= 0) {
        qname = bam_get_qname(read);
        if (mateseqs.count(qname) == 0 && mateseqs.count(qname + "_1") == 0 &&
            mateseqs.count(qname + "_2") == 0) continue; // mateseq not present

        std::string mate_read = mateseqs[qname];
        if (is_poly_ACGT(read) || is_poly_ACGT(mate_read.c_str(), mate_read.length())) continue;

        int pos = rev ? read->core.pos : bam_endpos(read);
        anchor_t a(rev ? 'L' : 'R', contig_id, pos, pos, rev ? is_left_clipped(read): is_right_clipped(read));
        cluster_t* c = new cluster_t(a, a, DISC_TYPES.DC, 1);
        c->id = clusters.size();
        clusters.push_back(c);
        clusters_map.insert(std::make_pair(c->a1.start, c));
        clusters_map.insert(std::make_pair(c->a1.end, c));
        reads.push_back(bam_dup1(read));
    }

    sam_itr_destroy(iter);
    bam_destroy1(read);

    // union-find datastructure
    int n_reads = reads.size();
    int* parents = new int[n_reads], * sizes = new int[n_reads];
    for (int i = 0; i < n_reads; i++) {
        parents[i] = i;
        sizes[i] = 1;
    }

    std::vector<int> max_dists;
    for (int i = 0; i < 10; i++) max_dists.push_back(i);
    for (int i = 10; i < 100; i+=10) max_dists.push_back(i);
    for (int i = 100; i < config.max_is; i+=100) max_dists.push_back(i);
    max_dists.push_back(config.max_is);

    for (int max_dist : max_dists) {
        std::priority_queue<cc_distance_t> pq;
        for (cluster_t* c1 : clusters) {
            if (c1->dead) continue;
            auto end = clusters_map.upper_bound(c1->a1.end+max_dist);
            for (auto map_it = clusters_map.lower_bound(c1->a1.start); map_it != end; map_it++) {
                cluster_t* c2 = map_it->second;
                if (c1 != c2 && cluster_t::can_merge(c1, c2, config) &&
                    (c1->a1.start <= c2->a1.start)) {
                    pq.push(cc_distance_t(cluster_t::distance(c1, c2), c1, c2));
                }
            }
        }

        while (!pq.empty()) {
            cc_distance_t ccd = pq.top();
            pq.pop();

            if (ccd.c1->dead || ccd.c2->dead) continue;

            cluster_t* new_cluster = cluster_t::merge(ccd.c1, ccd.c2);
            new_cluster->id = ccd.c1->id;
            merge(parents, sizes, ccd.c1->id, ccd.c2->id);
            clusters.push_back(new_cluster);


            ccd.c1->dead = true;
            remove_cluster_from_mm(clusters_map, ccd.c1);
            ccd.c2->dead = true;
            remove_cluster_from_mm(clusters_map, ccd.c2);

            auto end = clusters_map.upper_bound(new_cluster->a1.end + max_dist);
            for (auto map_it = clusters_map.lower_bound(new_cluster->a1.start - max_dist);
                 map_it != end; map_it++) {
                if (cluster_t::can_merge(new_cluster, map_it->second, config)) {
                    pq.push(cc_distance_t(cluster_t::distance(new_cluster, map_it->second), new_cluster,
                                          map_it->second));
                }
            }
            clusters_map.insert(std::make_pair(new_cluster->a1.start, new_cluster));
            clusters_map.insert(std::make_pair(new_cluster->a1.end, new_cluster));
        }
    }

    // for each set of reads, make a cluster-vector
    std::vector<std::vector<bam1_t*> > read_clusters;
    for (int i = 0; i < reads.size(); i++) {
        read_clusters.push_back(std::vector<bam1_t*>());
    }
    for (int i = 0; i < reads.size(); i++) {
        read_clusters[find(parents, i)].push_back(reads[i]);
    }

    // remap clusters
    std::vector<bam1_t*> kept_reads;
    for (int i = 0; i < read_clusters.size(); i++) {
        if (read_clusters[i].size() > 1) {
            int mq = 0;
            for (bam1_t* read : read_clusters[i]) { // check average MQ
                mq += read->core.qual;
            }
            remap_cluster(read_clusters[i], contig_id, header, rev, mateseqs, aligner, aligner_to_base, dc_remapped_file);
            for (bam1_t* read : read_clusters[i]) {
                kept_reads.push_back(read);
            }
        }
    }

    // write reads
    sort(kept_reads.begin(), kept_reads.end(), [](bam1_t *r1, bam1_t *r2) { return r1->core.pos < r2->core.pos; });
    for (bam1_t* r : kept_reads) {
        int ok = sam_write1(dc_remapped_file, header, r);
        if (ok < 0) throw "Unable to write to " + std::string(dc_remapped_file->fn);
    }

    for (cluster_t* c : clusters) delete c;

    delete[] parents; delete[] sizes;
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(dc_remapped_file);

    // destroy reads
    for (bam1_t* r : reads) {
        bam_destroy1(r);
    }

    dc_remapped_file = sam_open(dc_remapped_fname.c_str(), "r");

    code = sam_index_build(dc_remapped_file->fn, 0);
    if (code != 0) {
        throw "Cannot index " + std::string(dc_remapped_file->fn);
    }

    sam_close(dc_remapped_file);
    sam_close(dc_file);
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

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (contig_id = 1; contig_id < contig_id2name.size(); contig_id++) {
        std::future<void> future = thread_pool.push(remap, contig_id, true);
        futures.push_back(std::move(future));
        future = thread_pool.push(remap, contig_id, false);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    // split DC files
    const char* directions = "RL";
    for (contig_id = 1; contig_id < contig_id2name.size(); contig_id++) {
        for (int i = 0; i < 2; i++) {
            char dir = directions[i];
            std::vector<samFile *> r_writers(contig_id2name.size()+1, NULL), l_writers(contig_id2name.size()+1, NULL);

            char fname_str[1000];
            sprintf(fname_str, "%s/%c%d-DC.remap.bam", workspace.c_str(), dir, contig_id);
            samFile *reader = sam_open(fname_str, "r");

            std::cout << "Splitting " << fname_str << std::endl;

            hts_idx_t *idx = sam_index_load(reader, fname_str);
            if (idx == NULL) {
                throw "Unable to open index for " + std::string(fname_str);
            }

            bam_hdr_t *header = sam_hdr_read(reader);
            if (header == NULL) {
                throw "Unable to open header for " + std::string(fname_str);
            }

            std::string contig = contig_id2name[contig_id];
            hts_itr_t *iter = sam_itr_querys(idx, header, contig.c_str());
            bam1_t *read = bam_init1();

            while (sam_itr_next(reader, iter, read) >= 0) {
                int m_contig_id = contig_name2id[std::string(header->target_name[read->core.mtid])];
                char w_name[1000];
                if (bam_is_mrev(read)) {
                    if (l_writers[m_contig_id] == NULL) {
                        sprintf(w_name, "%s/%c%d-L%d-DC.bam", workspace.c_str(), dir, contig_id, m_contig_id);
                        l_writers[m_contig_id] = open_writer(std::string(w_name), header);
                    }
                    int ok = sam_write1(l_writers[m_contig_id], header, read);
                    if (ok < 0) throw "Could not write to " + std::string(w_name);
                } else {
                    if (r_writers[m_contig_id] == NULL) {
                        sprintf(w_name, "%s/%c%d-R%d-DC.bam", workspace.c_str(), dir, contig_id, m_contig_id);
                        r_writers[m_contig_id] = open_writer(std::string(w_name), header);
                    }
                    int ok = sam_write1(r_writers[m_contig_id], header, read);
                    if (ok < 0) throw "Could not write to " + std::string(w_name);
                }
            }

            for (int j = 1; j <= contig_id2name.size(); j++) {
                if (l_writers[j] != NULL) {
                    std::string writer_name = l_writers[j]->fn;
                    sam_close(l_writers[j]);
                    l_writers[j] = sam_open(writer_name.c_str(), "r");
                    int code = sam_index_build(l_writers[j]->fn, 0);
                    if (code != 0) {
                        throw "Cannot index " + std::string(l_writers[j]->fn);
                    }
                    sam_close(l_writers[j]);
                }
                if (r_writers[j] != NULL) {
                    std::string writer_name = r_writers[j]->fn;
                    sam_close(r_writers[j]);
                    r_writers[j] = sam_open(writer_name.c_str(), "r");
                    int code = sam_index_build(r_writers[j]->fn, 0);
                    if (code != 0) {
                        throw "Cannot index " + std::string(r_writers[j]->fn);
                    }
                    sam_close(r_writers[j]);
                }
            }
        }
    }
}