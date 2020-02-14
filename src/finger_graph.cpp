#include <zlib.h>
#include "gzstream.h"
#include "finger_graph.h"
#include "utils.h"

Lyndon::FingerGraph::FingerGraph(std::istream &in, int k, int limit, bool normalize, bool enriched_kfingers)
    : k(k), limit(limit), is_normalized(normalize), is_directed(!normalize), is_enriched(!enriched_kfingers) {
    std::string line; int c = 0;
    while (std::getline(in, line)) {
        auto split_line = split(line, '|');
        if (split_line.size() < 2) { continue; }

        auto read_id_offset = split(split_line[0], ' ');
        auto read_id = read_id_offset[0];
        auto offset = std::stoi(read_id_offset[1]);
        auto factors = split(split_line[1], ' ');
        auto fp = facts2fingerprint(factors);

        int idx = 0;
        while (idx + k < fp.size()) {
            auto kfL = slice(fp, idx, idx + k);
            auto key_seqL = enriched_kfingers ? Lyndon::get_key_factor(factors, idx, idx + k, this->is_normalized) : "";

            auto kfR = slice(fp, idx + 1, idx + k + 1);
            auto key_seqR = enriched_kfingers ? Lyndon::get_key_factor(factors, idx + 1, idx + k + 1, this->is_normalized) : "";

            if (sum(kfL) < limit || sum(kfR) < limit) {
                offset += fp[idx];
                idx++;
                continue;
            }

            Node* nL = add_node(
                    kfL, key_seqL, read_id, offset, join(factors, "", idx, idx + k)
            );
            Node* nR = add_node(
                    kfR, key_seqR, read_id, offset + kfL[0], join(factors, "", idx + 1, idx + k + 1)
            );
            add_edge(nL, nR);

            offset += fp[idx];
            idx++;
        }

    }
}

Lyndon::FingerGraph::FingerGraph() : k(0), limit(0), is_normalized(false), is_directed(true) { }

Lyndon::Node* Lyndon::FingerGraph::add_node(const Lyndon::k_finger &kf, const std::string &key_seq,
                                   const std::string &read_id, int offset, const std::string &seq) {
    auto key = make_key(kf, key_seq);
    if (this->nodes.find(key) == this->nodes.end()) {
        Node* n = make_node(kf, key_seq, read_id, offset, seq);
        this->nodes[key] = n;
        return n;
    } else {
        auto occ = Lyndon::Occurrence { read_id, offset };
        this->nodes[key]->occs.insert(occ);
        return this->nodes[key];
    }
}

void Lyndon::FingerGraph::add_edge(Lyndon::Node* n1, Lyndon::Node* n2) {
    n1->adj_list.insert(n2);
    if (!this->is_directed) {
        n2->adj_list.insert(n1);
    }
}

Lyndon::NodeKey Lyndon::FingerGraph::make_key(const Lyndon::k_finger &kf, const std::string &key_seq) const {
    if (this->is_normalized) {
        return Lyndon::NodeKey { normalize(kf), normalize(key_seq) };
    }
    return Lyndon::NodeKey { kf, key_seq };
}

Lyndon::Node* Lyndon::FingerGraph::make_node(const Lyndon::k_finger &kf, const std::string &key_seq,
                                             std::set<Lyndon::Occurrence> &occs) const {
    return new Lyndon::Node {
        make_key(kf, key_seq),
        occs,
        std::set<Node*>()
    };
}

Lyndon::Node* Lyndon::FingerGraph::make_node(const Lyndon::k_finger &kf, const std::string &key_seq,
                                             const std::string &read_id, int offset, const std::string &seq) const {

    auto node = new Lyndon::Node {
            make_key(kf, key_seq),
            std::set<Lyndon::Occurrence>(),
            std::set<Node*>()
    };
    auto occ = Lyndon::Occurrence { read_id, offset };
    node->occs.insert(occ);

    return node;
}

std::string Lyndon::get_key_factor(const Lyndon::factorization &factors, int begin, int end, bool normalize) {
    if (end - begin > 3) {
        begin += 1;
        end -= 1;
    }

    int idx = begin;
    int max_idx = idx, max_length = factors[idx].length();
    while (idx < end) {
        if (factors[idx].length() > max_length) {
            max_idx = idx;
            max_length = factors[idx].length();
        }

        idx++;
    }

    auto longest = factors[max_idx];
    if (normalize) {
        auto rev_comp = reverse_complement(longest);
        if (rev_comp < longest) {
            longest = rev_comp;
        }
    }

    if (longest.length() > 20) {
        longest = longest.substr(0, 10) + longest.substr(longest.length() - 10, 10);
    }
    return longest;
}

Lyndon::k_finger Lyndon::normalize(const Lyndon::k_finger &kf) {
    if (kf.size() == 0) {
        return kf;
    }

    int left = 0, right = kf.size() -1;
    while (left < right) {
        if (kf[left] < kf[right]) {
            return kf;
        }
        if (kf[left] > kf[right]) {
            auto result = std::vector<int>(kf.size());
            for (int i = 0; i < kf.size(); i++) {
                result[i] = kf[kf.size() - i - 1];
            }
            return result;
        }

        left++;
        right--;
    }

    return kf;
}

std::string Lyndon::normalize(const std::string &seq) {
    if (seq.length() == 0) {
        return seq;
    }

    int left = 0, right = seq.length() - 1;
    while (left < right) {
        if (seq[left] < complement(seq[right])) {
            return seq;
        }
        if (seq[left] > complement(seq[right])) {
            return reverse_complement(seq);
        }

        left++;
        right--;
    }

    return seq;
}

bool Lyndon::operator<(const Lyndon::Occurrence &x, const Lyndon::Occurrence &y) {
    if (x.r_id == y.r_id) {
        return x.offset < y.offset;
    }
    return x.r_id < y.r_id;
}

bool Lyndon::operator==(const Lyndon::Occurrence &x, const Lyndon::Occurrence &y) {
    return x.r_id == y.r_id && x.offset == y.offset;
}

std::ostream &operator<<(std::ostream &out, Lyndon::Occurrence o) {
    out << o.to_string();
    return out;
}

std::string Lyndon::Occurrence::to_string() const {
    return "('" + this->r_id + "', " + std::to_string(this->offset) + ")";
}

bool Lyndon::operator<(const Lyndon::NodeKey &x, const Lyndon::NodeKey &y) {
    if (x.kf == y.kf) {
        return x.key_sequence < y.key_sequence;
    }
    return x.kf < y.kf;
}

bool Lyndon::operator==(const Lyndon::NodeKey &x, const Lyndon::NodeKey &y) {
    return x.kf == y.kf && x.key_sequence == y.key_sequence;
}

bool Lyndon::operator<(const Lyndon::Node &x, const Lyndon::Node &y) {
    return x.key < y.key;
}

bool Lyndon::operator==(const Lyndon::Node &x, const Lyndon::Node &y) {
    return x.key == y.key;
}

Lyndon::FingerGraph* Lyndon::FingerGraph::from_graph_file(const std::string &file_path) {
    //    TODO: implement
    auto* graph = new Lyndon::FingerGraph();

    igzstream in(file_path.c_str());
    std::string line;

    getline(in, line);
    auto header_v = split(line, '\t');
    graph->k = std::stoi(header_v[1].substr(2));
    graph->limit = std::stoi(header_v[3].substr(10));
    graph->is_normalized = header_v[4].substr(14) == "True";
    graph->is_directed = !graph->is_normalized;

    while (getline(in, line)) {
        if (line.find("ED") == 0) break;
        if (line.find("VT") != 0) {
            fprintf(stderr, "Error while reading graph. File corrupted.\n");
            abort();
        }

        // Parse key
        auto kf = std::vector<int>();
        auto key_seq = "";
        // Parse occs set
        auto occs = std::set<Lyndon::Occurrence>();
        // Parse occ
        Lyndon::Occurrence occ = { };
        occs.insert(occ);

        auto new_node = graph->make_node(kf, key_seq, occs);
        graph->nodes[new_node->key] = new_node;

        std::cout << line << std::endl;
    }

    return graph;
}

void Lyndon::FingerGraph::save() {
    std::vector<std::string> htv = { "HT", "k=" + std::to_string(this->k),
        "threshold=" + std::to_string(this->limit), "is_normalized=" + std::to_string(this->is_normalized),
        "is_enriched=" + std::to_string(this->is_enriched) + "\n"};
    auto ht = join(htv, "\t");
    printf("%s", ht.c_str());

    for (const auto &pair : this->nodes) {
        const Node* n = pair.second;
        if (n->occs.size() < 2) {
            continue;
        }
        printf("VT\t((%s), '%s')\t%s\n",
                 v2s(n->key.kf, ", ").c_str(),
                n->key.key_sequence.c_str(),
                set2str(n->occs).c_str());
    }

    for (const auto &pair : this->nodes) {
        const Node* n1 = pair.second;

        if (n1->occs.size() < 2) {
            continue;
        }

        auto kf1_str = v2s(n1->key.kf, ", ");
        for (const Node* n2 : n1->adj_list) {
            if (n2->occs.size() < 2) {
                continue;
            }

            auto kf2_str = v2s(n2->key.kf, ", ");
            printf("ED\t((%s), '%s')\t((%s), '%s')\n",
                     kf1_str.c_str(),
                     n1->key.key_sequence.c_str(),
                     kf2_str.c_str(),
                     n2->key.key_sequence.c_str());
            if (!this->is_directed) {
                 printf("ED\t((%s), '%s')\t((%s), '%s')\n",
                     kf2_str.c_str(),
                     n2->key.key_sequence.c_str(),
                     kf1_str.c_str(),
                     n1->key.key_sequence.c_str());
            }
        }
    }
}

Lyndon::FingerGraph::~FingerGraph() {
    for (auto &pair : this->nodes) {
        delete pair.second;
    }
}
