#include <iostream>
#include <fstream>
#include <unordered_set>
#include <cassert>

#include "config.h"
#include "cluster.h"

config_t config;

std::vector<std::string> contig_id2name;

int support(prediction_t& pred) {
    return pred.disc_pairs + std::min(pred.bp1.sc_reads, pred.bp2.sc_reads);
}

double ptn_score(const prediction_t& pred) {
//    return std::max(double(pred.disc_pairs+pred.bp1.sc_reads)/std::max(pred.bp1.spanning_reads,1),
//                    double(pred.disc_pairs+pred.bp2.sc_reads)/std::max(pred.bp2.spanning_reads,1));
    return double(pred.disc_pairs+pred.bp1.sc_reads)/std::max(pred.bp1.spanning_reads,1);
}

std::string formatted_print(prediction_t& pred) {
    char str[10000];
    if (pred.bp1.sc_reads > 0 && pred.sv_type != SV_TYPES.INS) {
        sprintf(str, "ID=%d BP1=%c:%s:%d BP2=%c:%s:%d %s DISC=%d SC=%d,%d SPANNING=%d,%d PVAL=%lf EST_SIZE=%d SHIFT-PVAL=%lf",
                pred.id,
                pred.bp1.dir, contig_id2name[pred.bp1.contig_id].c_str(), pred.bp1.pos(),
                pred.bp2.dir, contig_id2name[pred.bp2.contig_id].c_str(), pred.bp2.pos(),
                svt_to_str(pred.sv_type).c_str(), pred.disc_pairs, pred.bp1.sc_reads, pred.bp2.sc_reads,
                pred.bp1.spanning_reads, pred.bp2.spanning_reads, pred.pval, pred.get_size(), pred.shift_pval);
    } else {
        sprintf(str, "ID=%d BP1=%c:%s:%d BP2=%c:%s:%d %s DISC=%d SC=%d,%d SPANNING=%d,%d PVAL=%lf EST_SIZE=%d:%d SHIFT-PVAL=%lf",
                pred.id,
                pred.bp1.dir, contig_id2name[pred.bp1.contig_id].c_str(), pred.bp1.pos(),
                pred.bp2.dir, contig_id2name[pred.bp2.contig_id].c_str(), pred.bp2.pos(),
                svt_to_str(pred.sv_type).c_str(), pred.disc_pairs, pred.bp1.sc_reads, pred.bp2.sc_reads,
                pred.bp1.spanning_reads, pred.bp2.spanning_reads, pred.pval, pred.size - pred.conf_ival,
                pred.size + pred.conf_ival, pred.shift_pval);
    }
    return str;
}

std::unordered_set<int> in_reps_ids;

bool operator < (const prediction_t& p1, const prediction_t& p2) {
//    bool in_rep1 = in_reps_ids.count(p1.id);
//    bool in_rep2 = in_reps_ids.count(p2.id);
//    if (in_rep1 != in_rep2) return in_rep1 < in_rep2;
//    if (p1.bp1.sc_reads*p1.bp2.sc_reads != p2.bp1.sc_reads*p2.bp2.sc_reads) {
//        return p1.bp1.sc_reads*p1.bp2.sc_reads > p2.bp1.sc_reads*p2.bp2.sc_reads;
//    }
//    if (p1.disc_pairs != p2.disc_pairs) {
//        return p1.disc_pairs > p2.disc_pairs;
//    }
//    if (p1.pval != p2.pval) {
//        return p1.pval < p2.pval;
//    }
    return ptn_score(p1) > ptn_score(p2);
    return p1.id < p2.id;
}

int main(int argc, char* argv[]) {
    std::string workdir = argv[1];
    bool no_filter = argc > 2 && std::string(argv[2]) == "no-filter";
    double ptn_ratio;
    int minsize;
    if (!no_filter) {
        ptn_ratio = argc > 2 ? std::stod(argv[2]) : 0.33;

        config = parse_config(workdir + "/config.txt");
    }

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
    }

    std::vector<prediction_t> retained;

    std::ifstream predictions_fin(workdir + "/predictions.info");
    std::string line;
    while (predictions_fin >> line) {
        prediction_t pred(line);

        if (no_filter) {
            retained.push_back(pred);
            continue;
        }

        if (support(pred) >= config.avg_depth/5 && ptn_score(pred) > ptn_ratio) {
            retained.push_back(pred);
        }
    }

    std::sort(retained.begin(), retained.end());
    for (prediction_t pred : retained) {
        std::cout << formatted_print(pred) << " " << in_reps_ids.count(pred.id) << std::endl;
    }
}
