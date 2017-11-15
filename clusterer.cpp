#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cassert>
#include <htslib/sam.h>

#include "config.h"
#include "sam_utils.h"
#include "libs/cptl_stl.h"
#include "cluster.h"

std::string workdir;
std::mutex mtx;

const int MAX_COV_FACTOR = 20;

std::vector<std::string> contig_id2name;
std::unordered_map<std::string, int> contig_name2id;

config_t config;

std::ofstream predictions_writer;


void clusterize(int id, int contig_id, int contig2_id, disc_type_t dt, std::string& fname, bool open_bam) {
    mtx.lock();
    std::cout << "Clustering " << fname << std::endl;
    mtx.unlock();

    std::string txt_fname = workdir + "/workspace/" + fname + ".txt";
    std::string bam_fname = workdir + "/workspace/" + fname + ".bam";

    std::vector<cluster_t*> clusters;

    samFile* bam_file = NULL;
    if (open_bam) bam_file = sam_open(bam_fname.c_str(), "r");
    if (bam_file != NULL) {
        int code = sam_index_build(bam_fname.c_str(), 0);
        if (code != 0) {
            throw "Cannot index " + bam_fname;
        }

        hts_idx_t* idx = sam_index_load(bam_file, bam_fname.c_str());
        if (idx == NULL) {
            throw "Unable to open index for " + bam_fname;
        }

        bam_hdr_t* header = sam_hdr_read(bam_file);
        if (header == NULL) {
            throw "Unable to open header for " + bam_fname;
        }

        std::string contig = contig_id2name[contig_id];
        hts_itr_t* iter = sam_itr_querys(idx, header, contig.c_str());
        bam1_t* read = bam_init1();
        while (sam_itr_next(bam_file, iter, read) >= 0) {
            anchor_t a1 = anchor_t(bam_is_rev(read) ? 'L' : 'R', contig_id, read->core.pos, bam_endpos(read), 0);
            anchor_t a2 = is_mate_unmapped(read) ? a1 :  // if mate is unmapped then make two copies of a1
                          anchor_t(bam_is_mrev(read) ? 'L' : 'R', contig2_id, read->core.mpos, get_mate_endpos(read), 0);
            clusters.push_back(new cluster_t(a1, a2, dt, 1));
        }
    }

    std::string line;
    std::ifstream txt_fin(txt_fname);
    while (txt_fin >> line) {
        clusters.push_back(new cluster_t(line, dt));
    }

    if (clusters.empty()) return;

    std::multimap<int, cluster_t*> clusters_map;
    for (cluster_t* cluster : clusters) {
        clusters_map.insert(std::make_pair(cluster->a1.start, cluster));
    }

    // remove extremely high coverage regions
    std::vector<cluster_t*> to_be_removed;
    int max_reads = (int) ((2*config.max_is * config.avg_depth * MAX_COV_FACTOR)/config.read_len);
    for (int i = 0; i < clusters.size(); i++) {
        cluster_t* c = clusters[i];
        int clusterable = 0;
        auto end = clusters_map.upper_bound(c->a1.end+config.max_is);
        for (auto map_it = clusters_map.lower_bound(c->a1.start-config.max_is); map_it != end; map_it++) {
            if (abs(map_it->second->a2.pos()-c->a2.pos()) <= config.max_is) {
                clusterable++;
                if (clusterable > max_reads) break;
            }
        }
        if (clusterable > max_reads) {
            clusters[i] = NULL;
            to_be_removed.push_back(c);
        }
    }
    clusters.erase(std::remove(clusters.begin(), clusters.end(), (cluster_t*) NULL), clusters.end());
    for (cluster_t* c : to_be_removed) {
        auto eq_range = clusters_map.equal_range(c->a1.start);
        for (auto it = eq_range.first; it != eq_range.second; it++) {
            if (it->second == c) {
                clusters_map.erase(it);
                break;
            }
        }
    }
    to_be_removed.clear();

    std::priority_queue<cc_distance_t> pq;
    for (cluster_t* c1 : clusters) {
        auto end = clusters_map.upper_bound(c1->a1.end+config.max_is);
        for (auto map_it = clusters_map.lower_bound(c1->a1.start); map_it != end; map_it++) {
            cluster_t* c2 = map_it->second;
            if (c1 != c2 && cluster_t::can_merge(c1, c2, config) &&
                (c1->a1.start <= c2->a1.start)) {
                pq.push(cc_distance_t(cluster_t::distance(c1, c2), c1, c2));
            }
        }
    }

    while (!pq.empty()) {
        cc_distance_t ccd = pq.top();
        pq.pop();

        if (ccd.c1->dead | ccd.c2->dead) continue;

        cluster_t* new_cluster = cluster_t::merge(ccd.c1, ccd.c2);
        if (new_cluster->dt == DISC_TYPES.LI && // stop if they try to merge into an invalid cluster
            new_cluster->a2.pos()-new_cluster->a1.pos() < 50) {
            continue;
        }
        clusters.push_back(new_cluster);

        // if memory is a problem, remove elements instead of just marking as dead
        ccd.c1->dead = true;
        ccd.c2->dead = true;

        auto end = clusters_map.upper_bound(new_cluster->a1.end+config.max_is);
        for (auto map_it = clusters_map.lower_bound(new_cluster->a1.start-config.max_is); map_it != end; map_it++) {
            if (cluster_t::can_merge(new_cluster, map_it->second, config)) {
                pq.push(cc_distance_t(cluster_t::distance(new_cluster, map_it->second), new_cluster, map_it->second));
            }
        }
        clusters_map.insert(std::make_pair(new_cluster->a1.start, new_cluster));
    }

    mtx.lock();
    for (cluster_t* c : clusters) {
        if (!c->dead) {
            prediction_t prediction(c, dt);
            predictions_writer << prediction.to_str() << "\n";
        }
    }
    mtx.unlock();
}

int main(int argc, char* argv[]) {
    workdir = std::string(argv[1]);
    std::string workspace = workdir + "/workspace";

    config = parse_config(workdir + "/config.txt");

    // we explicitly store the contig_name2id to make sure the order is consistent among all execs
    std::ifstream contig_map_fin(workdir + "/contig_map");
    std::string contig_name; int contig_id;
    contig_id2name.push_back("");
    while (contig_map_fin >> contig_name >> contig_id) {
        contig_id2name.push_back(contig_name);
        contig_name2id[contig_name] = contig_id;
    }

    predictions_writer.open(workspace + "/predictions.raw");

    std::ifstream dc_files_fin(workspace + "/dc-files.txt");
    std::string fname;
    std::vector<std::string> dc_bams;
    std::unordered_set<std::string> dc_txts;
    while (dc_files_fin >> fname) {
        if (fname[fname.length()-1] == 'm') { // ends in .bam
            dc_bams.push_back(fname.substr(0, fname.length()-7));
        } else if (fname[fname.length()-1] == 't') { // ends in .txt
            dc_txts.insert(fname.substr(0, fname.length()-7));
        }
    }

    ctpl::thread_pool thread_pool(config.threads);

    std::vector<std::future<void> > futures;
    for (std::string& dc_fname : dc_bams) {
        std::future<void> future;

        int contig2_id;
        char c;
        sscanf(dc_fname.c_str(), "%c%d-%c%d", &c, &contig_id, &c, &contig2_id);

        future = thread_pool.push(clusterize, contig_id, contig2_id, DISC_TYPES.DC, dc_fname + "-DC", true);
        futures.push_back(std::move(future));

        dc_txts.erase(dc_fname);
    }
    for (auto& it : dc_txts) {
        std::future<void> future;

        int contig2_id;
        char c;
        sscanf(it.c_str(), "%c%d-%c%d", &c, &contig_id, &c, &contig2_id);

        future = thread_pool.push(clusterize, contig_id, contig2_id, DISC_TYPES.DC, it + "-DC", false);
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
        } catch (std::string& s) {
            std::cout << s << std::endl;
        }
    }

    predictions_writer.close();
}
