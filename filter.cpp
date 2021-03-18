#include <iostream>
#include <fstream>
#include <unordered_set>

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

char dir_to_strand(char dir) {
    return (dir == 'R' ? 'F' : 'R');
}

std::string formatted_print(prediction_t& pred) {
    char str[10000];
    sprintf(str, "ID=%d BP1=%c:%s:%d BP2=%c:%s:%d DISCORDANT=%d SPLIT_READS=%d,%d SCORE=%lf",
            pred.id,
            dir_to_strand(pred.bp1.dir), contig_id2name[pred.bp1.contig_id].c_str(), pred.bp1.pos(),
            dir_to_strand(pred.bp2.dir), contig_id2name[pred.bp2.contig_id].c_str(), pred.bp2.pos(),
            pred.disc_pairs, pred.bp1.sc_reads, pred.bp2.sc_reads, ptn_score(pred));
    return str;
}

bool operator < (const prediction_t& p1, const prediction_t& p2) {
//    return ptn_score(p1) > ptn_score(p2);
    return p1.id < p2.id;
}

int main(int argc, char* argv[]) {
    std::string workdir = argv[1];
    bool no_filter = argc > 2 && std::string(argv[2]) == "no-filter";
    bool use_sc = false;
    double ptn_ratio;
    int minsize;
    if (!no_filter) {
        ptn_ratio = argc > 2 ? std::stod(argv[2]) : 0.33;
        if (argc > 3 && std::string(argv[3]) == "use_sc") use_sc = true;

        config = parse_config(workdir + "/config.txt");
    }

    // we explicitly store contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
//    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
    }

    std::vector<prediction_t> retained;

    std::ifstream predictions_fin(workdir + "/predictions.info");
    std::string line;
    while (predictions_fin >> line) {
        prediction_t pred(line);
        if (!use_sc) pred.bp1.sc_reads = pred.bp2.sc_reads = 0;

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
        std::cout << formatted_print(pred) << std::endl;
    }
}
