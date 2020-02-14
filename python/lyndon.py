import gzip
from tqdm import tqdm
from ast import literal_eval
from factorization import cfl, icfl, cfl_icfl, d_cfl, d_icfl, d_cfl_icfl, algs_dict
from utils import get_remove_border, reverse_complement, borders_dict


class DBG:
    @staticmethod
    def kfingers(factors, lengths, k):
        for i in range(len(lengths) - k + 1):
            yield tuple(factors[i:i+k]), tuple(lengths[i:i+k]) # f1f2f3
            # yield "$".join(factors[i:i+k]), tuple(lengths[i:i+k]) # f1$f2$f3

    @staticmethod
    def from_reads_file(filename, **kwargs):
        with open(filename) as f: data = f.read().split('\n')[:-1]
        return DBG(reads=data, **kwargs)

    @staticmethod
    def normalize(kf1):
        for f1, f2 in zip(kf1, reversed(kf1)):
            if f1 < f2:
                return kf1
            if f1 > f2:
                return tuple(reversed(kf1))
        return kf1
    
    @staticmethod
    def get_key_factor(factors, normalize=False):
        if len(factors) > 3:
            longest = max(factors[1:-1], key=lambda factor: len(factor))
        else:
            longest = max(factors, key=lambda factor: len(factor))

        if normalize:
            reverse = reverse_complement(longest)
            if reverse < longest:
                longest = reverse
        
        if 10 < len(longest) < 20:
            longest = longest[:10] + longest[-(len(longest) - 10):]
        elif len(longest) > 20:
            longest = longest[:10] + longest[-10:]
        
        return longest
    
    @staticmethod
    def from_file(filename):
        g = DBG(_empty=True)
        f = gzip.open(filename)

        type_, *params = f.readline().decode('utf-8').strip().split('\t')
        assert(type_ == "HT")
        assert(len(params) == 5)

        k_param, k_val = params[0].split('=')
        assert(k_param == "k")
        g.k = int(k_val)

        alg_param, alg_val = params[1].split('=')
        assert(alg_param == "alg")
        g.alg = algs_dict[alg_val]

        limit_param, param_val = params[2].split('=')
        assert(limit_param == "limit")
        g.limit = int(param_val)

        is_normalized_param, is_normalized_val = params[3].split('=')
        assert(is_normalized_param == "is_normalized")
        g.is_normalized = is_normalized_val == "True"
        g.is_directed = not g.is_normalized

        border_param, border_val = params[4].split('=')
        assert(border_param == "border_function")
        g.remove_border = borders_dict[border_val]

        while f:
            line = f.readline()
            type_, *rest = line.decode('utf-8').strip().split('\t')
            if type_ == "ED":
                break
            if type_ != "VT":
                raise Exception("File corrupted")

            key_str, value_str = rest

            key = literal_eval(key_str)
            seqs = literal_eval(value_str)
            node = DBG.Node(key[0], key[1])
            node.seqs = seqs

            g.nodes[key] = node
        
        while f:
            line = f.readline()
            if line == b'':
                break

            type_, key1_str, key2_str = line.decode('utf-8').strip().split('\t')

            if type_ != "ED":
                raise Exception("File corrupted")

            key1 = literal_eval(key1_str)
            key2 = literal_eval(key2_str)

            node1 = g.nodes[key1]
            node2 = g.nodes[key2]

            g.add_edge(node1, node2, directed=g.is_directed)
        
        return g
            
    
    class Node:
        def __init__(self, kf, longest_factor, seq=None, read_id=None, offset=None):
            self.key = (kf, longest_factor)
            self.kf = kf
            self.seqs = {}
            if seq is not None:
                self.seqs[(read_id, offset)] = seq

        def __hash__(self):
            return hash(self.key)
            # return hash(self.kf)

        def __str__(self):
            return "(" + " ".join(map(str, self.kf)) + "; Longest Factor: " + self.key[1] + ")"

    def __init__(self, reads=None, k=5, limit=1, remove_border=None, alg=None, normalize=False, _empty=False, **kwargs):
        if _empty:
            self.G = {}
            self.nodes = {}
            return
        if not alg:
            alg = cfl_icfl
            print("Using default cfl_icfl.")
        self.k = k
        self.alg = alg
        self.limit = limit
        self.is_normalized = normalize
        self.is_directed = not self.is_normalized
        self.remove_border = remove_border
        self.G = {}     # Adjacency list
        self.nodes = {} # k-fingers -> node

        k = k + 1
        for s in tqdm(reads):
            offset = 0
            read_id, read = s.split()
            factors = alg(read, **kwargs)
            lengths = [len(i) for i in factors]

            if remove_border is not None:
                l, lengths, r = remove_border(lengths)
                factors = factors[len(l):-len(r)]
                offset += sum(l)
            kfingers = list(self.kfingers(factors, lengths, k))

            for kf_factors, kf in kfingers:
                nodeL, nodeR = None, None
                kfL, kfR = kf[:-1], kf[1:]
                factorsL, factorsR = kf_factors[:-1], kf_factors[1:]

                if sum(kfL) < self.limit or sum(kfR) < self.limit:
                    offset += kf[0]
                    continue

                # f1f2f3
                seqL = ''.join(factorsL)
                seqR = ''.join(factorsR)
                # f1$f2$f3
                # seqL = '$'.join(factorsL)
                # seqR = '$'.join(factorsR)

                key_seqL = self.get_key_factor(factorsL, normalize=self.is_normalized)
                key_seqR = self.get_key_factor(factorsR, normalize=self.is_normalized)

                if normalize:
                    kfL = self.normalize(kfL)
                    kfR = self.normalize(kfR)

                nodeL = self.add_node(read_id, kfL, key_seqL, seqL, offset)
                nodeR = self.add_node(read_id, kfR, key_seqR, seqR, offset + kf[0])
                self.add_edge(nodeL, nodeR)

                offset += kf[0]


    def add_node(self, read_id, kf, longest_factor, seq, offset):
        key = (kf, longest_factor)
        if key in self.nodes:
            node = self.nodes[key]
            node.seqs[(read_id, offset)] = seq
        else:
            self.nodes[key] = self.Node(kf, longest_factor, seq, read_id, offset)
            node = self.nodes[key]
        return node

    def add_edge(self, nodeL, nodeR, directed=True):
        self.G.setdefault(nodeL, []).append(nodeR)
        if not self.is_directed:
            self.G.setdefault(nodeR, []).append(nodeL)

    def makeSet(self, x):
        return [x]

    def findSet(self, x, sets, set_member_lookup):
        return set_member_lookup[x]

    def union(self, setU, setV, sets, set_member_lookup):
        if sets[setU] is not None:
            if sets[setV] is not None:
                sets[setU].extend(sets[setV])

                for k in sets[setV]:
                    set_member_lookup[k] = setU
                sets[setV] = None

    def findCC(self):
        sets = []
        set_member_lookup = {}

        for index, key in enumerate(self.nodes):
            set_member_lookup[self.nodes[key]] = index
            sets.append(self.makeSet(self.nodes[key]))
        for left in tqdm(self.G):
            for right in self.G[left]:
                setU = self.findSet(left, sets, set_member_lookup)
                setV = self.findSet(right, sets, set_member_lookup)

                if setV != setU:
                    self.union(setU, setV, sets, set_member_lookup)

        return filter(None, sets)
