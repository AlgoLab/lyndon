#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <stdexcept>
#include <sys/types.h>
#include <random>
#include <algorithm>
#include <vector>
#include <utility>

#include "BooPHF.h"



std::vector<std::string>
split(const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string item;

    while (std::getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}


std::pair<std::vector<int>, std::pair<int, int>>
parse_line(const std::string &line) {
    auto splitted = split(line, ' ');
    std::vector<int> tokens;

    for (const std::string &str : splitted) {
        tokens.push_back(stoi(str));
    }

    auto id = tokens.back(); tokens.pop_back();
    auto pos = tokens.back(); tokens.pop_back();
    auto id_pos = std::make_pair(id, pos);

    return std::make_pair(tokens, id_pos);
}


std::pair<std::vector<std::vector<int>>, std::vector<std::pair<int, int>> >
parse_file(const std::string &path) {
    std::ifstream infile(path);
    std::string line;

    std::pair<std::vector<std::vector<int>>, std::vector<std::pair<int, int>> > result;

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        auto pair = parse_line(line);

        result.first.push_back(pair.first);
        result.second.push_back(pair.second);
    }

    return result;
};


void print_pairs(std::vector<std::vector<int>> facts, std::vector<std::pair<int, int>> pos_lens) {
    std::cout << "Size: " << facts.size() << std::endl << std::endl;

    for (auto i = 0; i < facts.size(); i++) {
        std::cout << facts[i].size() << "-ple: ";
        for (auto v : facts[i]) std::cout << v << " ";
        std::cout << std::endl << "Pair: " << pos_lens[i].first << pos_lens[i].second << std::endl << std::endl;
    }
}


class tuple_hasher
{
public:
//	uint64_t operator() (my_struct key, uint64_t seed=0) const
//	{
//		uint64_t hash[2];
//		MurmurHash3_x64_128(&key, sizeof(key), seed, hash);
//		return hash[0];
//	}

    uint64_t operator() (std::vector<int> key, uint64_t seed=0) const
    {
        uint64_t hash = 0;

        for (auto v : key) {
            hash ^= v >> 13;
        }

        hash ^= seed;

        return hash;
    }
};

typedef boomphf::mphf<std::vector<int>, tuple_hasher> boophf_t;



int main (int argc, char *argv[]) {
    auto pairs = parse_file(argv[1]);

    auto facts = pairs.first;
    auto pos_lens = pairs.second;

    if (argc > 2)
        print_pairs(facts, pos_lens);

    auto gammaFactor = 1.0;
    auto nthreads = 1;
    auto data_iterator = boomphf::range(facts.begin(), facts.end());
    auto* mpf = new boophf_t(
            facts.size(), data_iterator, nthreads, gammaFactor, false
    );

    std::vector<std::pair<int, int>> nodes_hash(facts.size());
    for (auto i = 0; i < facts.size(); i++) {
        auto idx = mpf->lookup(facts[i]);
        nodes_hash[idx] = pos_lens[i];
    }


}
