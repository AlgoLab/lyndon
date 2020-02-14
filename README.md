# Lyndon Finger Graph

## Build

```bash
g++ -std=c++11 -O3 -g ./src/*.cpp ./src/gzstream.C -I./src -o finger-graph -lz
```



## Run

### Compute Factorizations

```bash
usage: lfg compute-factorizations [-h] -a
                                  {cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb}
                                  [-b {remove-three,up-to-ten,twenty-most}]
                                  [-o O] [-n N]
                                  fasta

positional arguments:
  fasta                 fasta file (or folder if with -n)

optional arguments:
  -h, --help            show this help message and exit
  -a {cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb}
                        factorization algorithm
  -b {remove-three,up-to-ten,twenty-most}
                        strategy to remove borders
  -o O                  output file (or folder if with -n)
  -n N                  number of processes
```



#### Example

```bash
./lfg compute-factorizations -a cfl_icfl_comb reads.fa -o factorizations.txt
```



### Finger Graph

```bash
usage: lfg build [-h] [-k K] [-l LIMIT] [--no-norm] [--no-enriched]
                 factorizations

positional arguments:
  factorizations  factorizations file

optional arguments:
  -h, --help      show this help message and exit
  -k K            k-finger dimension [default 5]
  -l LIMIT        minimum length for a k-finger [default 30]
  --no-norm       do not normalize k-fingers
  --no-enriched   do not enrich k-fingers
```



#### Example

```bash
./lfg build -k 5 -l 30 factors.txt
./lfg build factors.txt | gzip > finger-graph.txt.gz
```
