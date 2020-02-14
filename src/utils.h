#ifndef LYNDON_UTILS_H
#define LYNDON_UTILS_H

#include <vector>
#include <string>
#include <sstream>
#include <stdarg.h>
#include "finger_graph.h"

void log(const char * format, ...);

std::tuple<int, int> remove_three(const Lyndon::fingerprint &fingerprint);
std::tuple<int, int> up_to_ten(const Lyndon::fingerprint &fingerprint);

std::string v2s(const std::vector<int> &v, const std::string &sep=" ");
std::string map2str(const std::map<Lyndon::Occurrence, std::string> &m);
std::string set2str(const std::set<Lyndon::Occurrence> &occs);
std::vector<std::string> split(const std::string &s, char delim);
std::string join(const std::vector<std::string> &v, std::string s, int start = 0, int end = -1);
std::vector<int> vecstr2vecint(const std::vector<std::string> &vs);

Lyndon::fingerprint facts2fingerprint(const std::vector<std::string> &facts);

std::string &ltrim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
std::string &rtrim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
std::string &trim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
int sum(const std::vector<int> &v, int start = 0, int end = -1);

std::map<Lyndon::read_id, Lyndon::factorization> load_factorizations(const std::string& path);
std::map<Lyndon::read_id, Lyndon::fingerprint> load_fingerprints(const std::string& path);

std::string reverse_complement(const std::string &seq);
inline char complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        case 'N':
            return 'N';
        default:
            abort();
    }
}

template <class T> std::vector<T> slice(const std::vector<T> &v, int start = 0, int end = -1) {
    {
        if (end == -1 || end > v.size()) {
            end = v.size();
        }
        if (start < 0) {
            start = 0;
        }

        if (end < start) {
            return std::vector<T>();
        }

        auto result = std::vector<T>(end - start);
        for (int i = start; i < end; i++) {
            result[i - start] = v[i];
        }

        return result;
    }
}

#endif
