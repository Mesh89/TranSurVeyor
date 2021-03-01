#include <fstream>
#include <iostream>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <unistd.h>
KSEQ_INIT(int, read)

#include "config.h"
#include "libs/cptl_stl.h"
#include "sam_utils.h"
#include "sw_utils.h"

config_t config;
std::mutex mtx;
std::string workdir;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, std::pair<char*, size_t> > chrs;

std::ofstream fa_out;

struct consensus_t {
    bool left_clipped;
    int contig_id, start, end;
    std::string consensus;
    int sc_reads;

    consensus_t(bool left_clipped, int contig_id, int start, int end, const std::string& consensus, int sc_reads)
            : left_clipped(left_clipped), contig_id(contig_id), start(start), end(end), consensus(consensus),
              sc_reads(sc_reads) {}

    char dir() { return left_clipped ? 'L' : 'R'; }

    std::string name() {
        std::stringstream ss;
        ss << contig_id << "_" << start << "_" << end << "_" << dir() << "_" << sc_reads;
        return ss.str();
    }

    int clipped_pos() const { return left_clipped ? start : end; }
};
std::vector<std::vector<consensus_t*> > consensus_vectors;

auto clip_start = [](bam1_t* r, bool left_clipped) {return left_clipped ? get_unclipped_start(r) : bam_endpos(r);};
auto clip_end = [](bam1_t* r, bool left_clipped) {return left_clipped ? r->core.pos : get_unclipped_end(r);};


std::string build_clip_consensus2(std::vector<bam1_t*>& clipped, bool left_clipped) {

    sort(clipped.begin(), clipped.end(), [&left_clipped](bam1_t* r1, bam1_t* r2) {
        return clip_start(r1, left_clipped) < clip_start(r2, left_clipped);});

    int first_clip_start = clip_start(clipped[0], left_clipped);
    hts_pos_t last_clip_end = 0;
    for (bam1_t* r : clipped) {
        last_clip_end = std::max(last_clip_end, clip_end(r, left_clipped));
    }
    int clip_len = last_clip_end - clip_start(clipped[0], left_clipped);

    if (clip_len < MIN_CLIP_CONSENSUS_LEN) return "";

    char consensus[MAX_READ_SUPPORTED];

    std::vector<int> clips_start, clips_end;
    std::vector<const uint8_t*> read_seqs;
    for (bam1_t* r : clipped) {
        clips_start.push_back(clip_start(r, left_clipped)-first_clip_start);
        clips_end.push_back(clip_end(r, left_clipped)-first_clip_start);
        read_seqs.push_back(bam_get_seq(r));
    }

    int s = 0;
    for (int i = 0; i < clip_len; i++) {
        while (s < clipped.size() && clips_end[s] < i) s++;

        int a = 0, c = 0, g = 0, t = 0;
        for (int j = s; j < clipped.size() && clips_start[j] <= i; j++) {
            if (clips_end[j] <= i) continue;

            bam1_t* r = clipped[j];
            int clip_start_in_read = left_clipped ? 0 : r->core.l_qseq - (get_unclipped_end(r)-bam_endpos(r));
            uint8_t nucl = bam_seqi(read_seqs[j], clip_start_in_read + i-clips_start[j]);
            if (nucl == 1) a++;
            else if (nucl == 2) c++;
            else if (nucl == 4) g++;
            else if (nucl == 8) t++;
        }

        if (a+c+g+t == 0) {
            return "";
        } else if (a >= c && a >= g && a >= t) {
            consensus[i] = 'A';
        } else if (c >= a && c >= g && c >= t) {
            consensus[i] = 'C';
        } else if (g >= a && g >= c && g >= t) {
            consensus[i] = 'G';
        } else if (t >= a && t >= c && t >= g) {
            consensus[i] = 'T';
        }
    }
    consensus[clip_len] = '\0';

    return std::string(consensus);
}

bool validate_clip(const consensus_t* consensus, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
    std::pair<char *, size_t> temp = chrs[contig_id2name[consensus->contig_id]];
    char* chr = temp.first;
    size_t chr_len = temp.second;

    int region_len = consensus->consensus.length() * 1.1;
    int region_start = consensus->left_clipped ? consensus->start-region_len : consensus->end;
    if (region_start < 0) region_start = 0;
    if (region_start+region_len >= chr_len) region_len = chr_len - region_start - 1;

    if (region_len <= 0) { // clipped at the end of a contig
        return true;
    }

    for (int i = region_start; i <= region_start+region_len; i++) {
        if (chr[i] == 'N') return false;
    }

    StripedSmithWaterman::Alignment alignemnt;
    int mask_len = consensus->consensus.length()/2;
    if (mask_len < 15) mask_len = 15; // suppress the warning
    try {
        aligner.Align(consensus->consensus.data(), chr+region_start, region_len, filter, &alignemnt, mask_len);
    } catch (std::bad_alloc e) {
        std::cout << "Sizes: " << consensus->consensus.size() << " " << region_len << std::endl;
    }
    if (alignemnt.query_end-alignemnt.query_begin >= consensus->consensus.length()*0.9
        && get_matches(alignemnt) >= consensus->consensus.length()*0.8
        && get_gaps(alignemnt) < alignemnt.query_end-alignemnt.query_begin*0.15) {
        return false; // clip can be realigned well next to the anchor
    } else {
        return true;
    }
}

consensus_t* build_clip_consensus(int contig_id, std::vector<bam1_t*>& clipped,
                                                         bool left_clipped, StripedSmithWaterman::Aligner& aligner,
                                                         StripedSmithWaterman::Filter& filter) {

    std::string consensus = build_clip_consensus2(clipped, left_clipped);
    if (consensus == "") return NULL;

    int first_clip_start = clip_start(clipped[0], left_clipped);
    std::vector<bam1_t*> accepted_clips;
    for (bam1_t* r : clipped) {
        const uint8_t* read_seq = bam_get_seq(r);
        int mm = 0;

        int clip_start_in_consensus = clip_start(r, left_clipped)-first_clip_start;
        int clip_end_in_consensus = clip_end(r, left_clipped)-first_clip_start;
        int clip_start_in_read = left_clipped ? 0 : r->core.l_qseq - (get_unclipped_end(r)-bam_endpos(r));

        for (int i = clip_start_in_consensus, j = clip_start_in_read; i < clip_end_in_consensus; i++, j++) {
            if (consensus[i] != get_base(read_seq, j)) {
                mm++;
            }
        }

        if (mm < std::ceil(MAX_SEQ_ERROR * consensus.length())) {
            accepted_clips.push_back(r);
        }
    }

    clipped.swap(accepted_clips);
    if (clipped.size() <= 1) {
        clipped.clear();
        return NULL;
    } else {
        hts_pos_t start = INT32_MAX, end = 0;
        for (bam1_t* r : clipped) {
            start = std::min(start, r->core.pos);
            end = std::max(end, bam_endpos(r));
        }
        consensus = build_clip_consensus2(clipped, left_clipped);
        return new consensus_t(left_clipped, contig_id, start, end, consensus, clipped.size());
    }
}

void write_clips_ids(std::ofstream& out, std::vector<bam1_t*>& clips) {
    for (bam1_t* clip : clips) {
        out << bam_get_qname(clip) << "\n";
    }
}

void build_clip_consensuses(int id, int contig_id) {
    mtx.lock();
    std::cout << "Building consensus for " << contig_id << "-CLIP.bam (" << contig_id2name[contig_id] << ")" << std::endl;
    mtx.unlock();

    StripedSmithWaterman::Aligner aligner(2,2,4,1,true);
    StripedSmithWaterman::Filter filter;

    std::vector<consensus_t*>& clip_consensuses = consensus_vectors[contig_id-1];

    std::string clipped_ids_fname = workdir + "/workspace/" + std::to_string(contig_id)+"-CLIPPED";
    std::ofstream clipped_ids_file(clipped_ids_fname);

    std::string bam_fname = workdir + "/workspace/" + std::to_string(contig_id)+"-CLIP.bam";
    samFile* bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file == NULL) {
        throw "Unable to open " + bam_fname;
    }

    hts_idx_t* idx = sam_index_load(bam_file, bam_fname.c_str());
    if (idx == NULL) {
        int code = sam_index_build(bam_fname.c_str(), 0);
        if (code != 0) {
            throw "Unable to open BAM index.";
        } else {
            idx = sam_index_load(bam_file, bam_fname.c_str());
        }
    }

    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        throw "Unable to open BAM header.";
    }

    hts_itr_t* iter = sam_itr_querys(idx, header, contig_id2name[contig_id].c_str());
    bam1_t* read = bam_init1();

    // divide soft-clipped reads into left-clipped and right-clipped
    std::vector<bam1_t*> lc_reads, rc_reads;
    while (sam_itr_next(bam_file, iter, read) >= 0) {
        if (is_left_clipped(read)) {
            lc_reads.push_back(bam_dup1(read));
        }
        if (is_right_clipped(read)) {
//            if (get_unclipped_end(read)-bam_endpos(read) >= MIN_CLIP_CONSENSUS_LEN) {
                rc_reads.push_back(bam_dup1(read));
//            }
        }
    }

    bam_destroy1(read);

    hts_itr_destroy(iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);

    sam_close(bam_file);

    auto lc_same_cluster = [](bam1_t* r1, bam1_t* r2) {return abs(r1->core.pos-r2->core.pos) <= config.max_sc_dist;};
    auto rc_same_cluster = [](bam1_t* r1, bam1_t* r2) {return abs(bam_endpos(r1)-bam_endpos(r2)) <= config.max_sc_dist;};

    std::vector<bam1_t*> curr_candidate_cluster;
    if (!lc_reads.empty()) {
        for (bam1_t* lc_read : lc_reads) {
            if (!curr_candidate_cluster.empty() &&
                !lc_same_cluster(curr_candidate_cluster[0], lc_read)) { // candidate cluster complete
                consensus_t* consensus = build_clip_consensus(contig_id, curr_candidate_cluster, true, aligner, filter);
                write_clips_ids(clipped_ids_file, curr_candidate_cluster);
                curr_candidate_cluster.clear();

                if (consensus != NULL) {
                    clip_consensuses.push_back(consensus);
                }
            }
            curr_candidate_cluster.push_back(lc_read);
        }
        // process last cluster
        consensus_t* consensus = build_clip_consensus(contig_id, curr_candidate_cluster, true, aligner, filter);
        write_clips_ids(clipped_ids_file, curr_candidate_cluster);
        curr_candidate_cluster.clear();
        if (consensus != NULL) {
            clip_consensuses.push_back(consensus);
        }
    }

    sort(rc_reads.begin(), rc_reads.end(), [](bam1_t* r1, bam1_t* r2) {return bam_endpos(r1) < bam_endpos(r2);});
    if (!rc_reads.empty()) {
        for (bam1_t* rc_read : rc_reads) {
            if (!curr_candidate_cluster.empty() &&
                !rc_same_cluster(curr_candidate_cluster[0], rc_read)) { // candidate cluster complete
                consensus_t* consensus = build_clip_consensus(contig_id, curr_candidate_cluster, false, aligner, filter);
                write_clips_ids(clipped_ids_file, curr_candidate_cluster);
                curr_candidate_cluster.clear();

                if (consensus != NULL) {
                    clip_consensuses.push_back(consensus);
                }
            }
            curr_candidate_cluster.push_back(rc_read);
        }
        // process last cluster
        consensus_t* consensus = build_clip_consensus(contig_id, curr_candidate_cluster, false, aligner, filter);
        write_clips_ids(clipped_ids_file, curr_candidate_cluster);
        curr_candidate_cluster.clear();
        if (consensus != NULL) {
            clip_consensuses.push_back(consensus);
        }
    }

    for (bam1_t* lc_read : lc_reads) {
        bam_destroy1(lc_read);
    }
    for (bam1_t* rc_read : rc_reads) {
        bam_destroy1(rc_read);
    }

    // remove non-validated clips
    clip_consensuses.erase(std::remove_if(clip_consensuses.begin(), clip_consensuses.end(),
                                          [&aligner, &filter](const consensus_t* c) {return !validate_clip(c, aligner, filter);}),
                           clip_consensuses.end());

    // sort by clipping position
    sort(clip_consensuses.begin(), clip_consensuses.end(), [](consensus_t* c1, consensus_t* c2) {
        return c1->clipped_pos() < c2->clipped_pos();});

    clipped_ids_file.close();
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

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
    }

    ctpl::thread_pool thread_pool(config.threads);

    for (contig_id = 1; contig_id < contig_id2name.size(); contig_id++) {
        consensus_vectors.push_back(std::vector<consensus_t *>());
    }

    std::vector<std::future<void> > futures;
    for (contig_id = 1; contig_id < contig_id2name.size(); contig_id++) {
        std::future<void> future = thread_pool.push(build_clip_consensuses, contig_id);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    fa_out.open(workspace + "/CLIPS.fa");
    for (std::vector<consensus_t*>& cv : consensus_vectors) {
        for (consensus_t* consensus : cv) {
            fa_out << ">" << consensus->name() << "\n" << consensus->consensus << "\n";
            delete consensus;
        }
    }
    fa_out.close();

    for (auto& chr : chrs) {
        delete[] chr.second.first;
    }
}
