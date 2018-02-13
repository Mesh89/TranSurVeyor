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
std::vector<std::ofstream*> mate_seqs_writers_by_tid; // indexed by tid, not contig_id

std::mutex mtx;

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


    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file == NULL) {
        throw "Unable to open BAM file.";
    }

    hts_idx_t* idx = sam_index_load(bam_file, bam_fname.c_str());
    if (idx == NULL) {
        throw "Unable to open BAM index.";
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    int contig_id = contig_name2id[contig];
    char region[1000];
    sprintf(region, "%s:%d-%d", contig.c_str(), 1, target_len);

    mtx.lock();
    std::cout << "Categorizing " << region << std::endl;
    mtx.unlock();

    hts_itr_t* iter = sam_itr_querys(idx, header, region);
    bam1_t* read = bam_init1();

    std::unordered_map<disc_type_t, samFile*> writers;

    samFile* clip_writer = get_writer(std::to_string(contig_id) + "-CLIP.bam", header);

    samFile* rdc_writer = get_writer("R" + std::to_string(contig_id) + "-DC.noremap.bam", header);
    samFile* ldc_writer = get_writer("L" + std::to_string(contig_id) + "-DC.noremap.bam", header);

    std::deque<bam1_t*> two_way_buffer, forward_buffer;

    int i = 0;
    while (sam_itr_next(bam_file, iter, read) >= 0 && i < MAX_BUFFER_SIZE-1) {
        if (is_valid(read, false)) {
            bam1_t* read2 = bam_dup1(read);
            add_to_queue(two_way_buffer, read2, 2*MAX_BUFFER_SIZE);
            add_to_queue(forward_buffer, read2, MAX_BUFFER_SIZE);
            i++;
        }
    }

    while (!forward_buffer.empty()) {
        while (sam_itr_next(bam_file, iter, read) >= 0) {
            if (is_valid(read, false)) {
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
        int64_t mq = get_mq(read);
        if (read->core.qual < MIN_MAPQ && mq < MIN_MAPQ) continue;

        // clipped read
        if (read->core.qual >= MIN_MAPQ && is_clipped(read) && check_SNP(read, two_way_buffer, config.avg_depth)) {
            int ok = sam_write1(clip_writer, header, read);
            if (ok < 0) throw "Failed to write to " + std::string(clip_writer->fn);
        }

        // we accept one of the mates having mapq 0 only if they are on different chromosomes
        if ((read->core.qual < MIN_MAPQ || mq < MIN_MAPQ)
            && !is_dc_pair(read) && !is_mate_unmapped(read)) {
            continue;
        }

        if (is_dc_pair(read)) {
            if (read->core.qual >= MIN_DC_MAPQ || mq >= MIN_DC_MAPQ) {
                if (read->core.qual >= mq && check_SNP(read, two_way_buffer, config.avg_depth)) { // stable end
                    if ((bam_is_rev(read) && !is_right_clipped(read)) || (!bam_is_rev(read) && !is_left_clipped(read))) {
                        int ok = sam_write1(bam_is_rev(read) ? ldc_writer : rdc_writer, header, read);
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
                    mtx.lock();
                    if (bam_is_rev(read)) {
                        rc(read_seq_chr);
                    }
                    std::string qname = bam_get_qname(read);
                    if (is_samechr(read)) {
                        if (read->core.isize > 0) qname += "_1";
                        else qname += "_2";
                    }
                    *mate_seqs_writers_by_tid[read->core.mtid] << qname << " ";
                    *mate_seqs_writers_by_tid[read->core.mtid] << read_seq_chr << "\n";
                    mtx.unlock();
                }
            }
        }
    }

    for (bam1_t* r : two_way_buffer) {
        bam_destroy1(r);
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);

    sam_close(clip_writer);

    sam_close(rdc_writer);
    sam_close(ldc_writer);

    sam_close(bam_file);
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

    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file == NULL) {
        std::cerr << "Unable to open BAM file." << std::endl;
        return -1;
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    tid_to_contig_id.resize(header->n_targets);
    for (int i = 0; i < header->n_targets; i++) {
        int contig_id = contig_name2id[std::string(header->target_name[i])];
        std::string fname = std::to_string(contig_id) + "-MATESEQS";
        mate_seqs_writers_by_tid.push_back(new std::ofstream(workdir + "/workspace/" + fname));
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
            std::cout << s << std::endl;
        }
    }

    for (int i = 0; i < header->n_targets; i++) {
        mate_seqs_writers_by_tid[i]->close();
        delete mate_seqs_writers_by_tid[i];
    }

    bam_hdr_destroy(header);
    sam_close(bam_file);
}
