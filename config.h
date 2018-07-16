#ifndef SURVEYOR_CONFIG_H
#define SURVEYOR_CONFIG_H

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

int MIN_MAPQ = 1;
int MIN_CLIP_LEN = 5;
int MIN_CLIP_CONSENSUS_LEN = 15;
int MIN_DC_MAPQ = 20;
double MAX_SEQ_ERROR = 0.04;

int MAX_READ_SUPPORTED = 10000;

struct config_t {
    int threads;
    int avg_depth;
    int max_is, min_is;
    int max_sc_dist, max_tra_size;
};


config_t parse_config(std::string file) {
    std::unordered_map<std::string, std::string> config_params;
    std::ifstream fin(file);
    std::string name, value;
    while (fin >> name >> value) {
        config_params[name] = value;
    }
    fin.close();

    config_t config;
    config.threads = stoi(config_params["threads"]);
    config.avg_depth = stoi(config_params["avg_depth"]);
    config.min_is = stoi(config_params["min_is"]);
    config.max_is = stoi(config_params["max_is"]);
    config.max_sc_dist = stoi(config_params["max_sc_dist"]);
    config.max_tra_size = stoi(config_params["max_tra_size"]);
    return config;
};

#endif //SURVEYOR_CONFIG_H
