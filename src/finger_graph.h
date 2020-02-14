#ifndef LYNDONHASH_H
#define LYNDONHASH_H

#include <unordered_map>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <map>


namespace Lyndon {
    typedef std::vector<int> k_finger;
    typedef std::vector<int> fingerprint;
    typedef std::vector<std::string> factorization;
    typedef std::string read_id;

    struct Occurrence {
        read_id r_id;
        int offset;

        std::string to_string() const;
    };
    bool operator<(const Occurrence &x, const Occurrence &y);
    bool operator==(const Occurrence &x, const Occurrence &y);
    std::ostream &operator<<(std::ostream &out, Occurrence o);

    struct NodeKey {
        k_finger kf;
        std::string key_sequence;
    };
    bool operator<(const NodeKey &x, const NodeKey &y);
    bool operator==(const NodeKey &x, const NodeKey &y);

    struct KeyHasher
    {
        std::size_t operator()(const NodeKey& k) const
        {
            using std::size_t;
            using std::hash;
            using std::string;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:
            auto result = hash<string>()(k.key_sequence);
            for (auto i : k.kf) {
                result ^= hash<int>()(i);
            }
            return result;
        }
    };

    struct Node {
        NodeKey key;
        std::set<Occurrence> occs;
        std::set<Node*> adj_list;
    };
    bool operator<(const Node &x, const Node &y);
    bool operator==(const Node &x, const Node &y);

    class FingerGraph {
    public:
        FingerGraph();
        FingerGraph(std::istream &factors, int k, int limit, bool normalize, bool enriched_kfingers);
        ~FingerGraph();

        int k;
        int limit;
        bool is_normalized;
        bool is_directed;
        bool is_enriched;

        std::unordered_map<NodeKey, Node*, KeyHasher> nodes; // (kf, seq) -> Node

        void save();
        static FingerGraph* from_graph_file(const std::string &file_path);

    private:
        NodeKey make_key(const k_finger &kf, const std::string &key_seq) const;
        Node* make_node(const k_finger &kf, const std::string &key_seq, std::set<Occurrence> &occs) const;
        Node* make_node(const k_finger &kf, const std::string &key_seq, const std::string &read_id, int offset, const std::string &seq) const;

        Node* add_node(const Lyndon::k_finger &kf, const std::string &key_seq, const std::string &read_id, int offset, const std::string &seq);
        void add_edge(Node* n1, Node* n2);
    };

    std::string get_key_factor(const factorization &factors, int begin, int end, bool normalize);
    k_finger normalize(const k_finger &kf);
    std::string normalize(const std::string &seq);
}

#endif //LYNDONHASH_H
