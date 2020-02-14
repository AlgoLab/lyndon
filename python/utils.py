import random, sys
from datetime import datetime


def log(message, files=[sys.stderr]):
    dt = datetime.now()
    for file in files:
        print(f"[{dt:%H:%M}] - {message}", file=file)


def read_file(filename):
    with open(filename) as f:
        return f.read().split('\n')[:-1]

def read_dict_bytes(filename, value_int=False):
    result = dict()

    with open(filename) as f:
        for line in f:
            read_id, *value = line.rstrip().split(' ')
            read_id_b = read_id.encode('utf-8')
            if value_int:
                value = [int(x) for x in value]
            result[read_id_b] = value

    return result

def random_sequence(length):
    ab = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(ab) for _ in range(length))

def sequence_factors_lengths(length, alg):
    s = random_sequence(length)
    factors = alg(s)
    lengths = [len(f) for f in factors]
    return s, factors, lengths

def overlap(begin1, end1, begin2, end2):
    return begin1 <= end2 and begin2 <= end1

def get_overlap(begin1, end1, begin2, end2):
    return min(end1, end2) - max(begin1, begin2) + 1

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[i] for i in reversed(seq))

def get_remove_border(border):
    def f(l):
        return l[:border], l[border:-border], l[-border:]
    return f

def sample_read_no_error(region, coverage=37.5,
        read_len=150, reverse_and_complement=True):
    reads = []
    n_read = int(len(region) * coverage / read_len)
    for _ in range(n_read):
        b = random.randint(0, len(region)-read_len)
        seq = region[b:b+read_len]
        flip = 0
        if reverse_and_complement:
            flip = random.randint(0,1)
            if flip: seq = reverse_complement(seq)
        reads.append(f"{b+1}_{b+150}_{flip} {seq}")
    return reads

def print_connected_components(ccs):
    for i, c in enumerate(ccs):
        print(f"Component {i}, len: {len(c)}")
        for n in c:
            print(f"{n.kf} -> [", end='')
            for read_id, offset in n.seqs:
                print(f"({read_id}, {offset}), ", end='')
            print("]")
        print("\n")


def write_connected_components(ccs, outfile):
    f = open(outfile, 'w')
    for i, c in enumerate(ccs):
        f.write(f"Component {i}, len: {len(c)}\n")
        for n in c:
            f.write(f"{n.kf} -> [")
            for read_id, offset in n.seqs:
                f.write(f"({read_id}, {offset}), ")
            f.write("]\n")
        f.write("\n\n")


def levenshtein(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def save_graph_to_file(g, out_file):
    out_file.write("\t".join(["HT", f"k={g.k}", f"alg={g.alg.__name__}", f"limit={g.limit}", f"is_normalized={g.is_normalized}", f"border_function={g.remove_border.__name__}\n"]))
    for node in g.nodes.values():
        key = node.key
        if len(node.seqs) < 2: continue

        out_file.write(f"VT\t({key[0]}, '{key[1]}')\t")
        out_file.write(str(node.seqs))
        out_file.write("\n")

    for node in g.nodes.values():
        key = node.key
        node_keys = list(node.seqs.keys())
        if len(node_keys) < 2: continue

        adj_list = g.G[node]
        for adj in adj_list:
            if len(adj.seqs) < 2: continue
            adj_key = adj.key
            out_file.write(f"ED\t({key[0]}, '{key[1]}')\t({adj_key[0]}, '{adj_key[1]}')\n")


def custom_border(lengths):
    l = lengths[:2]
    if (sum(l) >= 20):
        l = lengths[:1]
    r = lengths[-2:]
    if (sum(r) >= 20):
        r = lengths[-1:]
    m = lengths[len(l):-len(r)]

    return l, m, r

twenty_most = custom_border

def remove_three(lengths):
    if len(lengths) <= 6:
        return [], [], []
    return lengths[:3], lengths[3:-3], lengths[-3:]

def up_to_ten(lengths):
    if len(lengths) < 3:
        return [], [], []
    m = lengths.copy()
    l = [m.pop(0)]; r = [m.pop()]
    while m and sum(l) + m[0] <= 10:
        l.append(m.pop(0))
    while m and sum(r) + m[-1] <= 10:
        r.append(m.pop())
    return l, m, r

borders = [custom_border, remove_three, up_to_ten]
borders_dict = {
    f.__name__: f for f in borders
}