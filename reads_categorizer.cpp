#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include "sam_utils.h"
#include "config.h"
#include "libs/cptl_stl.h"

config_t config;
std::unordered_map<std::string, int> contig_name2id;
std::vector<int> tid_to_contig_id;
std::vector<std::string> mate_seqs_fnames; // indexed by tid, not contig_id
std::vector<std::vector<std::string> > mate_seqs;

std::mutex mtx;
std::mutex* mtx_contig;

const int MAX_BUFFER_SIZE = 100;


std::string workdir;

bam1_t* add_to_queue(std::deque<bam1_t*>& q, bam1_t* o, int size_limit) {
    bam1_t* t = NULL;
    while (q.size() >= size_limit) {
        t = q.front();
        q.pop_front();
    }
    q.push_back(o);
    return t;
}

samFile* get_writer(std::string name, bam_hdr_t* header) {
    samFile* writer = sam_open((workdir + "/workspace/" + name).c_str(), "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + (workdir + name);
    }
    return writer;
}

void categorize(int id, std::string contig, std::string bam_fname, int target_len) {

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());

    int contig_id = contig_name2id[contig];
    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), 1, target_len);

    mtx.lock();
    std::cout << "Categorizing " << region << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region);
    bam1_t* read = bam_init1();

    std::unordered_map<disc_type_t, samFile*> writers;

    samFile* clip_writer = get_writer(std::to_string(contig_id) + "-CLIP.bam", bam_file->header);

    samFile* rdc_writer = get_writer("R" + std::to_string(contig_id) + "-DC.noremap.bam", bam_file->header);
    samFile* ldc_writer = get_writer("L" + std::to_string(contig_id) + "-DC.noremap.bam", bam_file->header);

    std::deque<bam1_t*> two_way_buffer, forward_buffer;

    int i = 0;
    while (sam_itr_next(bam_file->file, iter, read) >= 0 && i < MAX_BUFFER_SIZE-1) {
        if (is_valid(read, false, true)) {
            bam1_t* read2 = bam_dup1(read);
            add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
            add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
            i++;
        }
    }

    while (!forward_buffer.empty()) {
        while (sam_itr_next(bam_file->file, iter, read) >= 0) {
            if (is_valid(read, false, true)) {
                bam1_t* read2 = bam_dup1(read);
                bam1_t* to_destroy = add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
                bam_destroy1(to_destroy);
                add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
                break;
            }
        }

        bam1_t* read = forward_buffer.front();
        forward_buffer.pop_front();

        bam_aux_get(read, "MQ");
        uint8_t mapq = is_primary(read) ? read->core.qual : 60;
        int64_t mq = get_mq(read);
        if (mapq < MIN_MAPQ && mq < MIN_MAPQ) continue;

        // clipped read
        if (mapq >= MIN_MAPQ && is_clipped(read) && check_SNP(read, two_way_buffer, config.avg_depth)) {
            int ok = sam_write1(clip_writer, bam_file->header, read);
            if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        }

        // we accept one of the mates having mapq 0 only if they are on different chromosomes
        if ((mapq < MIN_MAPQ || mq < MIN_MAPQ)
            && !is_dc_pair(read) && !is_mate_unmapped(read)) {
            continue;
        }

        if (is_dc_pair(read)) {
            if (mapq >= MIN_DC_MAPQ || mq >= MIN_DC_MAPQ) {
                if (mapq >= mq && check_SNP(read, two_way_buffer, config.avg_depth)) { // stable end
                    if ((bam_is_rev(read) && !is_right_clipped(read)) || (!bam_is_rev(read) && !is_left_clipped(read))) {
                        int ok = sam_write1(bam_is_rev(read) ? ldc_writer : rdc_writer, bam_file->header, read);
                        if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
                    }
                }
                if (mapq <= mq) { // save read seq for remapping
                    const uint8_t* read_seq = bam_get_seq(read);
                    char read_seq_chr[MAX_READ_SUPPORTED];
                    for (int i = 0; i < read->core.l_qseq; i++) {
                        read_seq_chr[i] = get_base(read_seq, i);
                    }
                    read_seq_chr[read->core.l_qseq] = '\0';
                    if (bam_is_rev(read)) {
                        rc(read_seq_chr);
                    }
                    std::string qname = bam_get_qname(read);
                    if (is_samechr(read)) {
                        if (read->core.isize > 0) qname += "_1";
                        else qname += "_2";
                    }
                    mtx_contig[read->core.mtid].lock();
                    mate_seqs[read->core.mtid].push_back(qname + " " + read_seq_chr);
                    mtx_contig[read->core.mtid].unlock();
                }
            }
        }
    }

    for (bam1_t* r : two_way_buffer) {
        bam_destroy1(r);
    }

    bam_destroy1(read);

    sam_close(clip_writer);

    sam_close(rdc_writer);
    sam_close(ldc_writer);

    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {
    std::string bam_fname = argv[1];
    workdir = std::string(argv[2]);
    std::string workspace = workdir + "/workspace";

    config = parse_config(workdir + "/config.txt");

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_name2id[contig_name] = contig_id;
    }

    ctpl::thread_pool thread_pool(config.threads);

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());

    bam_hdr_t* header = bam_file->header;
    tid_to_contig_id.resize(header->n_targets);
    mtx_contig = new std::mutex[header->n_targets + 1];
    mate_seqs.resize(header->n_targets + 1);
    for (int i = 0; i < header->n_targets; i++) {
        int contig_id = contig_name2id[std::string(header->target_name[i])];
        std::string fname = std::to_string(contig_id) + "-MATESEQS";
        mate_seqs_fnames.push_back(workdir + "/workspace/" + fname);
    }

    std::vector<std::future<void> > futures;
    for (int i = 0; i < header->n_targets; i++) {
        std::future<void> future = thread_pool.push(categorize, header->target_name[i], bam_fname, header->target_len[i]);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (char const* s) {
            std::cerr << s << std::endl;
        }
    }

    for (int i = 0; i < header->n_targets; i++) {
        std::ofstream mate_seqs_fout(mate_seqs_fnames[i]);
        for (std::string& mate_seq : mate_seqs[i]) {
            mate_seqs_fout << mate_seq << std::endl;
        }
        mate_seqs_fout.close();
    }

    close_samFile(bam_file);
}

