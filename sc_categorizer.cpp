#include <iostream>
#include <fstream>
#include <vector>
#include <htslib/sam.h>

#include "libs/cptl_stl.h"
#include "config.h"
#include "sam_utils.h"
#include "cluster.h"

config_t config;
std::mutex mtx;
std::string workdir;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;

bam_hdr_t* header;

void categorize(int id, int contig_id, std::string& clip_fname) {
    mtx.lock();
    std::cout << "Categorizing SC for " << contig_id << " (" << contig_id2name[contig_id] << ")" << std::endl;
    mtx.unlock();

    samFile* bam_file = sam_open(clip_fname.c_str(), "r");
    if (bam_file == NULL) {
        throw "Unable to open BAM file.";
    }

    hts_idx_t* idx = sam_index_load(bam_file, clip_fname.c_str());
    if (idx == NULL) {
        throw "Unable to open BAM index.";
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    std::vector<std::ofstream*> rr_dc_writers, rl_dc_writers, lr_dc_writers, ll_dc_writers;
    for (int i = 0; i <= contig_id2name.size(); i++) {
        rr_dc_writers.push_back(NULL);
        rl_dc_writers.push_back(NULL);
        lr_dc_writers.push_back(NULL);
        ll_dc_writers.push_back(NULL);
    }

    std::string contig = contig_id2name[contig_id];
    hts_itr_t* iter = sam_itr_querys(idx, header, contig.c_str());
    bam1_t* read = bam_init1();

    while (sam_itr_next(bam_file, iter, read) >= 0) {
        // clip must nearly all (>80%) aligned
        if (!is_valid(read, false) || read->core.qual < MIN_MAPQ || bam_endpos(read)-read->core.pos < read->core.l_qseq*0.8) continue;

        int anchor_contig_id, anchor_start, anchor_end; char anchor_dir; int anchor_sc_reads;
        char* seq_name = bam_get_qname(read);
        sscanf(seq_name, "%d_%d_%d_%c_%d", &anchor_contig_id, &anchor_start, &anchor_end, &anchor_dir, &anchor_sc_reads);

        auto opp_dir = [] (char dir) { return dir == 'L' ? 'R' : 'L'; };
        anchor_t a_anchor(anchor_dir, anchor_contig_id, anchor_start, anchor_end, anchor_sc_reads);
        anchor_t a_clip(bam_is_rev(read) ? anchor_dir : opp_dir(anchor_dir),
                        contig_id, read->core.pos, bam_endpos(read), 1);

        // TODO: add UM case

        std::ofstream* writer;
        bool anchor_first;
        disc_type_t dt;
        if (anchor_contig_id != contig_id || abs(read->core.pos-anchor_start) > 100000) {
            anchor_first = true;
            dt = DISC_TYPES.DC;

            // choose correct DC writer among RR, RL, LR, LL
            if (anchor_dir == 'L' && !bam_is_rev(read)) {
                if (ll_dc_writers[anchor_contig_id] == NULL) {
                    ll_dc_writers[anchor_contig_id] = new std::ofstream(
                            workdir + "/workspace/L" + std::to_string(anchor_contig_id) + "-L" + std::to_string(contig_id) + "-DC.txt"
                    );
                }
                writer = ll_dc_writers[anchor_contig_id];
            } else if (anchor_dir == 'L' && bam_is_rev(read)) {
                if (lr_dc_writers[anchor_contig_id] == NULL) {
                    lr_dc_writers[anchor_contig_id] = new std::ofstream(
                            workdir + "/workspace/L" + std::to_string(anchor_contig_id) + "-R" + std::to_string(contig_id) + "-DC.txt"
                    );
                }
                writer = lr_dc_writers[anchor_contig_id];
            } else if (anchor_dir == 'R' && bam_is_rev(read)) {
                if (rl_dc_writers[anchor_contig_id] == NULL) {
                    rl_dc_writers[anchor_contig_id] = new std::ofstream(
                            workdir + "/workspace/R" + std::to_string(anchor_contig_id) + "-L" + std::to_string(contig_id) + "-DC.txt"
                    );
                }
                writer = rl_dc_writers[anchor_contig_id];
            } else if (anchor_dir == 'R' && !bam_is_rev(read)) {
                if (rr_dc_writers[anchor_contig_id] == NULL) {
                    rr_dc_writers[anchor_contig_id] = new std::ofstream(
                            workdir + "/workspace/R" + std::to_string(anchor_contig_id) + "-R" + std::to_string(contig_id) + "-DC.txt"
                    );
                }
                writer = rr_dc_writers[anchor_contig_id];
            }
        } else continue;

        cluster_t cluster(anchor_first ? a_anchor : a_clip, anchor_first ? a_clip : a_anchor, dt, 0);
        *writer << cluster.to_str() << "\n";
    }

    for (int i = 1; i <= contig_id2name.size(); i++) {
        if (rr_dc_writers[i] != NULL) rr_dc_writers[i]->close();
        if (rl_dc_writers[i] != NULL) rl_dc_writers[i]->close();
        if (lr_dc_writers[i] != NULL) lr_dc_writers[i]->close();
        if (ll_dc_writers[i] != NULL) ll_dc_writers[i]->close();
    }
}

int main(int argc, char* argv[]) {

    workdir = std::string(argv[1]);
    std::string workspace = workdir + "/workspace";
    std::string clip_fname = workspace + "/CLIPS.sorted.bam";
    samFile* clip_file = sam_open(clip_fname.c_str(), "r");

    int code = sam_index_build(clip_fname.c_str(), 0);
    if (code != 0) {
        throw "Cannot read " + clip_fname;
    }

    header = sam_hdr_read(clip_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    config = parse_config(workdir + "/config.txt");

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (int contig_id = 1; contig_id <= header->n_targets; contig_id++) {
        std::future<void> future = thread_pool.push(categorize, contig_id, clip_fname);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
}

bam1_t* add_to_queue(std::deque<bam1_t*>& q, bam1_t* o, int size_limit) {
    bam1_t* t = NULL;
    while (q.size() >= size_limit) {
        t = q.front();
        q.pop_front();
    }
    q.push_back(o);
    return t;
}