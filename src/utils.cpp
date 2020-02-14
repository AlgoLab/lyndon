#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "utils.h"
#include <ctime>

void log(const char * format, ...) {
    time_t t = time(0);
    tm *ltm = localtime(&t);

    printf("[%d:%d] - ", ltm->tm_hour, ltm->tm_min);
    va_list args;
    va_start(args, format);
//    vfprintf(stderr, format, args);
    vprintf(format, args);
    va_end(args);
    flush(std::cout);
}

std::tuple<int, int> remove_three(const Lyndon::fingerprint &fingerprint) {
    auto size = fingerprint.size();
    if (size > 6) {
        return std::make_tuple(3, size - 3);
    } else {
        return std::make_tuple(0, 0);
    }
}

std::tuple<int, int> up_to_ten(const Lyndon::fingerprint &fingerprint) {
    auto size = fingerprint.size();
    if (size < 3) {
        return std::make_tuple(0, 0);
    }

    unsigned long left = 1, right = size - 2, sum = fingerprint[0];
    while (left < size && sum + fingerprint[left] <= 10) {
        sum += fingerprint[left];
        left++;
    }

    sum = fingerprint[size - 1];
    while (right >= 0 && sum + fingerprint[right] <= 10) {
        sum += fingerprint[right];
        right--;
    }

    return std::make_tuple(left, right + 1);

}

std::string v2s(const std::vector<int> &v, const std::string &sep)
{
    std::string res = "";
    for (std::size_t i = 0; i < v.size(); i++)
    {
        res += std::to_string(v[i]);
        if (i < v.size() - 1)
            res += sep;
    }
    return res;
}

std::string map2str(const std::map<Lyndon::Occurrence, std::string> &m) {
    std::string result = "{";

    for (auto it = m.begin(); it != m.end(); ) {
        result += it->first.to_string();
        result += ": ";
        result += "'" + it->second + "'";

        it++;

        if (it != m.end()) {
            result += ", ";
        }
    }

    result += "}";

    return result;
}

std::string set2str(const std::set<Lyndon::Occurrence> &occs) {
    std::string result = "{";

    for (auto it = occs.begin(); it != occs.end(); ) {
        result += (*it).to_string();
        it++;
        if (it != occs.end()) {
            result += ", ";
        }
    }
    result += "}";

    return result;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, delim))
    {
        result.push_back(item);
    }

    return result;
}

std::vector<int> vecstr2vecint(const std::vector<std::string> &vs)
{
    std::vector<int> result;
    for (auto s : vs)
    {
        try {
            result.push_back(stoi(s));
        } catch (...) {
            result.push_back(-2);
        }
    }
    return result;
}

Lyndon::fingerprint facts2fingerprint(const std::vector<std::string> &facts) {
    std::vector<int> result(facts.size());

    for (int i = 0; i < result.size(); i++) {
        result[i] = facts[i].length();
    }

    return result;
}

std::string &ltrim(std::string &str, const std::string &chars)
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}

std::string &rtrim(std::string &str, const std::string &chars)
{
    auto idx = str.find_last_not_of(chars);
    if (idx < str.length() - 1)
        str.erase(idx + 1);
    return str;
}

std::string &trim(std::string &str, const std::string &chars)
{
    if (str.length() > 0)
        return ltrim(rtrim(str, chars), chars);
    else
        return str;
}

int sum(const std::vector<int> &v, int start, int end) {
    if (end == -1 || end > v.size()) {
        end = v.size();
    }
    if (start < 0) {
        start = 0;
    }

    if (end < start || start >= v.size()) {
        return 0;
    }

    int sum = 0;
    for (int i = start; i < end; i++) {
        sum += v[i];
    }
    return sum;
}

std::string join(const std::vector<std::string> &v, std::string s, int start, int end) {
    if (v.empty()) {
        return "";
    }

    if (end == -1 || end > v.size()) {
        end = v.size();
    }
    if (start < 0) {
        start = 0;
    }

    if (end < start || start >= v.size()) {
        return "";
    }

    std::string result = v[start];
    for (int i = start + 1; i < end; i++) {
        result += s;
        result += v[i];
    }
    return result;
}

std::string reverse_complement(const std::string &seq) {
    auto result = std::string(seq.rbegin(), seq.rend());
    for (int i = 0; i < result.length(); i++) {
        result[i] = complement(result[i]);
    }
    return result;
}

std::map<Lyndon::read_id, Lyndon::factorization> load_factorizations(const std::string& path) {
    std::ifstream in(path);
    auto result = std::map<Lyndon::read_id, Lyndon::factorization>();

    std::string line;
    while (std::getline(in, line)) {
        auto read_id_length = line.find_first_of(' ');
        auto read_id = line.substr(0, read_id_length);

        auto idx = line.find_first_of('|') + 1;
        auto cur_start = idx;
        auto factors = std::vector<std::string>();
        while (idx < line.length()) {
            while (idx < line.length() && line[idx] != ' ') {
                idx += 1;
            }
            factors.push_back(line.substr(cur_start, idx - cur_start));
            idx += 1;
            cur_start = idx;
        }
        result[read_id] = factors;
    }

    return result;
}

std::map<Lyndon::read_id, Lyndon::fingerprint> load_fingerprints(const std::string& path) {
    std::ifstream in(path);
    auto result = std::map<Lyndon::read_id, Lyndon::fingerprint>();

    std::string line;
    while (std::getline(in, line)) {
        auto read_id_length = line.find_first_of(' ');
        auto read_id = line.substr(0, read_id_length);

        auto idx = line.find_first_of('|') + 1;
        auto cur_start = idx;
        auto fingerprint = std::vector<int>();
        while (idx < line.length()) {
            while (idx < line.length() && line[idx] != ' ') {
                idx += 1;
            }
            fingerprint.push_back(std::stoi(line.substr(cur_start, idx - cur_start)));
            idx += 1;
            cur_start = idx;
        }
        result[read_id] = fingerprint;
    }

    return result;
}
