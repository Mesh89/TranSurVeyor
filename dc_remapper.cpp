#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <unistd.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <cassert>
#include <htslib/hts.h>

KSEQ_INIT(int, read)

#include "sam_utils.h"
#include "config.h"
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
const double BASE_ACCEPTANCE_THRESHOLD = 0.9;

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
        if (alignment.sw_score >= r->core.l_qseq) {
            score += alignment.sw_score;
        }

        if (offsets != NULL) {
            if (alignment.sw_score < r->core.l_qseq) {
                offsets->push_back(SKIP_READ);
            } else {
                offsets->push_back(alignment.ref_begin);
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
                   StripedSmithWaterman::Aligner& aligner, samFile* dc_remapped_file) {
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

    compute_score(base_region, cluster, mateseqs, NULL, NULL, aligner, filter, is_rc);

    if (base_region.score >= best_region.score*BASE_ACCEPTANCE_THRESHOLD) {
        return;
    }

    // reorder by position - this way the DC files will not need to be sorted
    sort(cluster.begin(), cluster.end(), [] (bam1_t* r1, bam1_t* r2) {return r1->core.pos < r2->core.pos;});

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
        int ok = sam_write1(dc_remapped_file, header, r);
        if (ok < 0) throw "Unable to write to " + std::string(dc_remapped_file->fn);
    }
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

    StripedSmithWaterman::Aligner aligner(2,2,3,1);

    std::string contig = contig_id2name[contig_id];
    hts_itr_t* iter = sam_itr_querys(idx, header, contig.c_str());
    bam1_t* read = bam_init1();

    std::vector<bam1_t*> cluster;
    while (sam_itr_next(dc_file, iter, read) >= 0) {
        qname = bam_get_qname(read);
        if (mateseqs.count(qname) == 0 && mateseqs.count(qname + "_1") == 0 &&
            mateseqs.count(qname + "_2") == 0) continue; // mateseq not present

        if (is_poly_ACGT(read) || is_poly_ACGT(mateseqs[qname])) continue;

        if (!cluster.empty() && bam_endpos(read)-cluster[0]->core.pos > config.max_is) {
            remap_cluster(cluster, contig_id, header, rev, mateseqs, aligner, dc_remapped_file);
            cluster.clear();
        }
        cluster.push_back(bam_dup1(read));
    }
    remap_cluster(cluster, contig_id, header, rev, mateseqs, aligner, dc_remapped_file);

    sam_close(dc_remapped_file);
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
//        std::future<void> future = thread_pool.push(remap, contig_id, true);
//        futures.push_back(std::move(future));
        contig_id = 7;
        std::future<void> future = thread_pool.push(remap, contig_id, false);
        futures.push_back(std::move(future));
        break;
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    return 0;

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