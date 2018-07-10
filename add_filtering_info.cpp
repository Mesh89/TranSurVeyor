#include <iostream>
#include <fstream>
#include <queue>
#include <unordered_set>
#include <numeric>

#include "htslib/sam.h"
#include "config.h"
#include "cluster.h"
#include "libs/cptl_stl.h"

std::string workdir;
std::mutex mtx;

config_t config;

int MAX_READ_IS;

const int MAX_BUFFER_SIZE = 100;

std::unordered_map<std::string, int> contig_name2tid;
std::unordered_map<std::string, size_t> contig_name2len;
std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;

std::queue<open_samFile_t*> bam_pool;

open_samFile_t* get_bam_reader(std::string bam_fname) {
    if (!bam_pool.empty()) {
        open_samFile_t* o = bam_pool.front();
        bam_pool.pop();
        return o;
    }

    open_samFile_t* o = open_samFile(bam_fname.c_str());
    return o;
}

void release_bam_reader(open_samFile_t* reader) {
    bam_pool.push(reader);
}


int POP_SIZE = 10000;

bam1_t* add_to_queue(std::deque<bam1_t*>& q, bam1_t* o, int size_limit) {
    bam1_t* t = NULL;
    while (q.size() >= size_limit) {
        t = q.front();
        q.pop_front();
    }
    q.push_back(o);
    return t;
}


void make_contig_name2tid(std::string bam_fname) {
    open_samFile_t* open_sam = open_samFile(bam_fname.c_str());
    samFile* bam_file = open_sam->file;
    bam_hdr_t* header = open_sam->header;

    for (int i = 0; i < header->n_targets; i++) {
        contig_name2tid[std::string(header->target_name[i])] = i;
        contig_name2len[std::string(header->target_name[i])] = header->target_len[i];
    }
    close_samFile(open_sam);
}


void find_spanning(breakpoint_t& bp, std::string& bam_fname) {
    std::string contig = contig_id2name[bp.contig_id];

    mtx.lock();
    open_samFile_t* open_sam = get_bam_reader(bam_fname);
    mtx.unlock();

    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), bp.pos()-2*config.max_is, bp.pos());

    hts_itr_t* iter = sam_itr_querys(open_sam->idx, open_sam->header, region);
    bam1_t* read = bam_init1();
    while (sam_itr_next(open_sam->file, iter, read) >= 0) {
        if (!is_valid(read, false) || bam_is_rev(read)) continue;

        if (is_samechr(read) && !is_samestr(read) && !is_outward(read, config.min_is) &&
            read->core.isize >= config.min_is && read->core.isize <= config.max_is) {
            if ((bp.dir == 'R' && bp.pos() >= read->core.pos && bp.pos() < read->core.mpos) ||
                (bp.dir == 'L' && bp.pos() > bam_endpos(read) && bp.pos() <= get_mate_endpos(read))) {
                bp.spanning_reads++;
            }
        }
    }

    mtx.lock();
    release_bam_reader(open_sam);
    mtx.unlock();
};

void odds_ratio(int id, std::string bam_fname, prediction_t* pred) {
    find_spanning(pred->bp1, bam_fname);
    find_spanning(pred->bp2, bam_fname);
}

int main(int argc, char* argv[]) {
    std::string bam_fname = argv[1];
    workdir = argv[2];
    std::string workspace = workdir + "/workspace";

    config = parse_config(workdir + "/config.txt");

    std::vector<prediction_t*> preds;

    std::ifstream predictions_fin(workspace + "/predictions.raw");
    std::string line;
    while (predictions_fin >> line) {
        prediction_t* pred = new prediction_t(line);
        if (pred->disc_pairs+pred->bp1.sc_reads > 1 && pred->disc_pairs+pred->bp2.sc_reads > 1) {
            preds.push_back(pred);
        } else {
            delete pred;
        }
    }

    make_contig_name2tid(bam_fname);
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
    for (prediction_t* pred : preds) {
        std::future<void> future = thread_pool.push(odds_ratio, bam_fname, pred);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    std::ofstream predictions_filter_out(workdir + "/predictions.info");
    for (prediction_t* pred : preds) {
        predictions_filter_out << pred->to_str() << "\n";
        delete pred;
    }
    predictions_filter_out.close();

    while (!bam_pool.empty()) {
        close_samFile(bam_pool.front());
        bam_pool.pop();
    }
}
