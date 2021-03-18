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
std::string workdir;
const int MAX_BUFFER_SIZE = 100;

std::vector<std::vector<std::string> > mate_seqs;
std::mutex mtx;

std::mutex* mtx_contig;

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

void categorize(int id, int contig_id, std::string contig_name, std::string bam_fname) {

    mtx.lock();
    std::cout << "Categorizing " << contig_name << std::endl;
    mtx.unlock();

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());

    samFile* clip_writer = get_writer(std::to_string(contig_id) + "-CLIP.bam", bam_file->header);
    samFile* rdc_writer = get_writer("R" + std::to_string(contig_id) + "-DC.noremap.bam", bam_file->header);
    samFile* ldc_writer = get_writer("L" + std::to_string(contig_id) + "-DC.noremap.bam", bam_file->header);

    std::deque<bam1_t*> two_way_buffer, forward_buffer;

    int i = 0;
    bam1_t* read = bam_init1();
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, contig_name.c_str());
    while (sam_itr_next(bam_file->file, iter, read) >= 0 && i < MAX_BUFFER_SIZE-1) {
        if (is_valid(read)) {
            bam1_t* read2 = bam_dup1(read);
            add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
            add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
            i++;
        }
    }

    while (!forward_buffer.empty()) {
        while (sam_itr_next(bam_file->file, iter, read) >= 0) {
            if (is_valid(read)) {
                bam1_t* read2 = bam_dup1(read);
                bam1_t* to_destroy = add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
                bam_destroy1(to_destroy);
                add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
                break;
            }
        }

        bam1_t* read = forward_buffer.front();
        forward_buffer.pop_front();

        // clipped read
        if (read->core.qual >= config.min_stable_mapq && is_clipped(read, MIN_CLIP_LEN) && check_SNP(read, two_way_buffer, config.avg_depth)) {
            int ok = sam_write1(clip_writer, bam_file->header, read);
            if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        }

        int64_t mq = get_mq(read);
        if (is_dc_pair(read) && (read->core.qual >= config.min_stable_mapq || mq >= config.min_stable_mapq)) {
            if (read->core.qual >= mq && check_SNP(read, two_way_buffer, config.avg_depth)) { // stable end
                if ((bam_is_rev(read) && !is_right_clipped(read, MIN_CLIP_LEN)) ||
                    (!bam_is_rev(read) && !is_left_clipped(read, MIN_CLIP_LEN))) {
                    int ok = sam_write1(bam_is_rev(read) ? ldc_writer : rdc_writer, bam_file->header, read);
                    if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
                }
            }
            if (read->core.qual <= mq) { // save read seq for remapping
                const uint8_t* read_seq = bam_get_seq(read);
                char read_seq_chr[MAX_READ_SUPPORTED];
                for (int i = 0; i < read->core.l_qseq; i++) {
                    read_seq_chr[i] = get_base(read_seq, i);
                }
                read_seq_chr[read->core.l_qseq] = '\0';
                if (bam_is_rev(read)) {
                    rc(read_seq_chr);
                }
                std::string qname = bam_get_qname(read), read_seq_str = get_sequence(read);
                if (bam_is_rev(read)) {
                    rc(read_seq_str);
                }
                if (is_samechr(read)) {
                    if (read->core.isize > 0) qname += "_1";
                    else qname += "_2";
                }
                mtx_contig[read->core.mtid].lock();
                mate_seqs[read->core.mtid].push_back(qname + " " + read_seq_str);
                mtx_contig[read->core.mtid].unlock();
            }
        }
    }

    for (bam1_t* r : two_way_buffer) {
        bam_destroy1(r);
    }
    bam_destroy1(read);
    hts_itr_destroy(iter);

    sam_close(clip_writer);
    sam_close(rdc_writer);
    sam_close(ldc_writer);
    close_samFile(bam_file);
}

int main(int argc, char* argv[]) {
    std::string bam_fname = argv[1];
    workdir = std::string(argv[2]);
    std::string workspace = workdir + "/workspace/";

    config = parse_config(workdir + "/config.txt");

    open_samFile_t* bam_file = open_samFile(bam_fname.c_str());
    bam_hdr_t* header = bam_file->header;

    mtx_contig = new std::mutex[header->n_targets];
    mate_seqs.resize(header->n_targets);

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    for (int i = 0; i < header->n_targets; i++) {
        std::future<void> future = thread_pool.push(categorize, i, header->target_name[i], bam_fname);
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
        std::string fname = std::to_string(i) + "-MATESEQS";
        std::ofstream mate_seqs_fout(workspace + fname);
        for (std::string& mate_seq : mate_seqs[i]) {
            mate_seqs_fout << mate_seq << std::endl;
        }
        mate_seqs_fout.close();
    }

    close_samFile(bam_file);
}

