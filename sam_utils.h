#ifndef SURVEYOR_SAM_UTILS_H
#define SURVEYOR_SAM_UTILS_H

#include <iostream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstring>

#include "htslib/sam.h"

extern int MIN_MAPQ;
extern int MIN_CLIP_LEN;
extern double MAX_SEQ_ERROR;
extern int MAX_READ_SUPPORTED;

int get_mate_endpos(const bam1_t* r);

typedef uint8_t disc_type_t;
static struct disc_types_t {
    disc_type_t UM, DC, SS, OW, LI, SI;
    size_t n_types = 6;

    disc_types_t() : UM(1), DC(2), SS(3), OW(4), LI(5), SI(6) {}
} DISC_TYPES;

std::string dt_to_str(disc_type_t dt) {
    if (dt == DISC_TYPES.UM) return "UM";
    else if (dt == DISC_TYPES.DC) return "DC";
    else if (dt == DISC_TYPES.SS) return "SS";
    else if (dt == DISC_TYPES.OW) return "OW";
    else if (dt == DISC_TYPES.LI) return "LI";
    else if (dt == DISC_TYPES.SI) return "SI";
    else return "";
}

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}


bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}

bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid;
}

bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000;
}

bool is_samestr(bam1_t* r) {
    return (r->core.flag & BAM_FREVERSE) == (r->core.flag & BAM_FMREVERSE);
}

// assume pair is NOT same stranded, and this is leftmost read in pair
bool is_outward(bam1_t* r, int min_is) {
    if ((r->core.flag & BAM_FREVERSE) == 0) { // read if positive-strand
        int endpos = bam_endpos(r), m_endpos = get_mate_endpos(r);
        return (endpos > r->core.mpos && r->core.isize < min_is) || endpos > m_endpos;
    } else { // read is negative-strand
        return r->core.pos < r->core.mpos;
    }
}

// r must be to the left of mate
bool is_overlaps_with_mate(bam1_t* r) {
    return r->core.mpos < bam_endpos(r);
}

int64_t get_mq(bam1_t* r) {
    uint8_t* mq = bam_aux_get(r, "MQ");
    return mq == NULL ? 0 : bam_aux2i(mq);
}

bool is_valid(bam1_t* r) {
    return !is_unmapped(r) && is_primary(r);
}

bool is_left_clipped(bam1_t* r, int min_clip_len) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' && bam_cigar_oplen(cigar[0]) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' && bam_cigar_oplen(cigar[r->core.n_cigar-1]) >= min_clip_len;
}
bool is_clipped(bam1_t* r, int min_clip_len) {
    return is_left_clipped(r, min_clip_len) || is_right_clipped(r, min_clip_len);
}

int get_unclipped_start(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return r->core.pos - (bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]) : 0);
}

int get_unclipped_end(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_endpos(r) + (bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' ? bam_cigar_oplen(cigar[r->core.n_cigar-1]) : 0);
}

int get_endpoint(bam1_t* r) {
    return bam_is_rev(r) ? r->core.pos : bam_endpos(r);
}

std::string get_cigar_code(bam1_t* r) {
    const uint32_t* cigar = bam_get_cigar(r);
    std::stringstream ss;
    for (int i = 0; i < r->core.n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}

bool is_first_in_pair(bam1_t* r, disc_type_t dt) {
    if (dt == DISC_TYPES.UM) return true;
    else if (dt == DISC_TYPES.DC) {
        return r->core.qual >= get_mq(r);
    } else if (dt == DISC_TYPES.SS || dt == DISC_TYPES.SI || dt == DISC_TYPES.LI) {
        return r->core.pos < r->core.mpos;
    } else if (dt == DISC_TYPES.OW) {
        return r->core.flag & BAM_FREVERSE;
    } else return false;
}

int get_mate_endpos(const bam1_t* r) {
    uint8_t *mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0, pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}

int get_read_pos_given_ref_pos(bam1_t* r, int ref_pos) {
    if (r->core.pos > ref_pos || bam_endpos(r) < ref_pos) return -1;

    uint32_t* cigar = bam_get_cigar(r);

    int i = 0;
    int curr_pos_in_read = 0, curr_pos_in_ref = r->core.pos;
    if (bam_cigar_opchr(cigar[0]) == 'S') {
        curr_pos_in_read = bam_cigar_oplen(cigar[0]);
        i++;
    }

    for ( ; i < r->core.n_cigar; i++) {
        if (curr_pos_in_ref == ref_pos) return curr_pos_in_read;

        char opchr = bam_cigar_opchr(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        int move = std::min(ref_pos - curr_pos_in_ref, oplen);
        if (opchr == 'M') {
            curr_pos_in_read += move;
            curr_pos_in_ref += move;
        } else if (opchr == 'D') {
            curr_pos_in_ref += move;
        } else if (opchr == 'I') {
            curr_pos_in_read += oplen;
        }
    }

    return curr_pos_in_read;
}

bool check_SNP(bam1_t* read, std::deque<bam1_t*>& buffer, double avg_depth) {
    while (!buffer.empty() && bam_endpos(buffer.front()) < read->core.pos) {
        bam_destroy1(buffer.front());
        buffer.pop_front();
    }

    int* tot = new int[read->core.l_qseq];
    std::fill(tot, tot+read->core.l_qseq, 0);
    int* diff = new int[read->core.l_qseq];
    std::fill(diff, diff+read->core.l_qseq, 0);

    hts_pos_t endpos = bam_endpos(read);
    for (bam1_t* r : buffer) {
        if (r->core.pos > endpos) break;

        int us = std::max(read->core.pos, r->core.pos);
        int ue = std::min(endpos, bam_endpos(r));
        int i = get_read_pos_given_ref_pos(read, us);
        int ei = get_read_pos_given_ref_pos(read, ue);
        int j = get_read_pos_given_ref_pos(r, us);
        int ej = get_read_pos_given_ref_pos(r, ue);

        if (i < 0 || j < 0) continue;

        const uint8_t* read_seq = bam_get_seq(read);
        const uint8_t* r_seq = bam_get_seq(r);
        for (; i < ei && j < ej; i++, j++) {
            tot[i]++;

            uint8_t c1 = bam_seqi(read_seq, i);
            uint8_t c2 = bam_seqi(r_seq, j);

            if (c1 != c2) diff[i]++;
        }
    }

    int bonus = (int) std::ceil(MAX_SEQ_ERROR * (endpos-read->core.pos+1));
    for (int i = 0; i < read->core.l_qseq; i++) {
        if (tot[i] > 0 && tot[i]-diff[i] < std::min(avg_depth*0.4, tot[i]*0.4)) {
            bonus--;
        }

        if (bonus < 0) {
            delete[] tot;
            delete[] diff;
            return false;
        }
    }
    delete[] tot;
    delete[] diff;
    return true;
}

char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}

std::string get_sequence(bam1_t* r) {
    char seq[MAX_READ_SUPPORTED];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    return std::string(seq);
}

bool is_poly_ACGT(bam1_t* r) {
    int a = 0, c = 0, g = 0, t = 0;
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        char base = get_base(bam_seq, i);
        if (base == 'A') a++;
        else if (base == 'C') c++;
        else if (base == 'G') g++;
        else if (base == 'T') t++;
    }

    int maxc = std::max(std::max(a,c), std::max(g,t));
    return double(maxc)/r->core.l_qseq >= 0.8;
}

bool is_poly_ACGT(const char* seq, int len) {
    int a = 0, c = 0, g = 0, t = 0;
    for (int i = 0; i < len; i++) {
        char base = seq[i];
        if (base == 'A' || base == 'a') a++;
        else if (base == 'C' || base == 'c') c++;
        else if (base == 'G' || base == 'g') g++;
        else if (base == 'T' || base == 't') t++;
    }
    int maxc = std::max(std::max(a,c), std::max(g,t));
    return double(maxc)/len >= 0.8;
}

void rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

void rc(char* read) {
    int len = strlen(read);
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}


struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() {}

    open_samFile_t(samFile* file, bam_hdr_t* header, hts_idx_t* idx) : file(file), header(header), idx(idx) {}
};

open_samFile_t* open_samFile(const char* fname, bool index_file = false) {
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
    if (sam_file->idx == NULL) {
        throw "Unable to open index for " + std::string(fname);
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
    hts_idx_destroy(f->idx);
    bam_hdr_destroy(f->header);
    sam_close(f->file);
    delete f;
}


#endif //SURVEYOR_SAM_UTILS_H
