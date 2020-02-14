#!/usr/bin/env python3
import sys, fileinput, argparse, os, multiprocessing
from functools import partial
from os import path
from Bio import SeqIO
from tqdm import tqdm
from factorization import cfl, icfl, cfl_icfl, d_cfl, d_icfl, d_cfl_icfl
from utils import twenty_most, remove_three, up_to_ten

def compute_factorizations(file_path, output_path, alg):
    input_file = SeqIO.parse(file_path, 'fasta')
    if output_path is None:
        output_file = sys.stdout
    else:
        base_name = path.splitext(path.basename(file_path))[0]
        output_file = open(path.join(output_path, f"{base_name}.txt"), 'w')

    for r in input_file:
        rid, seq = r.id, str(r.seq)
        factors = alg(seq)
        output_file.write(rid + ' 0|')
        output_file.write(' '.join(factors) + '\n')

def compute_factorizations_with_border(file_path, output_path, alg, remove_border):
    input_file = SeqIO.parse(file_path, 'fasta')
    if output_path is None:
        output_file = sys.stdout
    else:
        base_name = path.splitext(path.basename(file_path))[0]
        output_file = open(path.join(output_path, f"{base_name}.txt"), 'w')

    for r in input_file:
        rid, seq = r.id, str(r.seq)
        factors = alg(seq)
        lengths = [len(x) for x in factors]
        l, _, r = remove_border(lengths)
        mid = factors[len(l):-len(r)]

        if len(mid) == 0:continue
        output_file.write(f"{rid} {sum(l)}|")
        output_file.write(' '.join(mid) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alg', dest='alg', action='store', default='cfl')
    parser.add_argument('-o', '--out', dest='output_path', action='store', default=None)
    parser.add_argument('-i', '--input', dest='input_path', action='store', default=None)
    parser.add_argument('--border', dest='border', action='store', default=None)
    parser.add_argument('-n', dest='n_processes', action='store', default=1, type=int)

    args = parser.parse_args()

    algs = {
        'cfl': cfl,
        'icfl': icfl,
        'cfl_icfl': cfl_icfl,
        'cfl_comb': d_cfl,
        'icfl_comb': d_icfl,
        'cfl_icfl_comb': d_cfl_icfl
    }

    borders = {
        'remove-three': remove_three,
        'up-to-ten': up_to_ten,
        'twenty-most': twenty_most
    }

    alg = algs[args.alg]

    input_path = args.input_path
    output_path = args.output_path
    input_path_is_dir = False
    if path.isdir(input_path):
        filenames = [f for f in os.listdir(input_path) if path.isfile(path.join(input_path, f))]
        input_path_is_dir = True
        if output_path is not None and not path.exists(output_path):
            os.makedirs(output_path)
    else:
        filenames = [input_path]
        dir_name = path.dirname(output_path)
        if not path.exists(dir_name):
            os.makedirs(dir_name)

    pool = multiprocessing.Pool(args.n_processes)
    if args.border is not None:
        border = borders[args.border]
        func = partial(compute_factorizations_with_border, alg=alg, remove_border=border, output_path=output_path)
    else:
        func = partial(compute_factorizations, alg=alg, output_path=output_path)

    for _ in pool.imap_unordered(func, [path.join(input_path, f) for f in filenames]):
        pass
