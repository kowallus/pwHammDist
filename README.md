# pwHammDist - pairwise Hamming distance algorithms
 
The repository contains a set of practical algorithms for finding pairs of 
similar sequences in a collection (the similarity measure is the Hamming distance).
The software works with sequences represented in unsigned 16-bit integers or as 
a set of basic DNA symbols. Benchmarks include also generation and matching of
synthetic binary data.

### Installation on Linux

The following steps create benchmark executables. 
On Linux PgRC build requires installed cmake version >= 3.4 
(check using ```cmake --version```):
```bash
git clone https://github.com/kowallus/pwHammDist.git
cd pwHammDist
mkdir build
cd build
cmake ..
make bench-binary
make bench-dataset
```

### Benchmarking datasets - basic usage 

```
bench-dataset [-a algorithmID] [-n] [-p] [-R] [-g] [-s] 
		[-r noOfRepeats] [-v] [-q] [-A] [-c] [-i] [-I] [-Q] datasetFileName k

algorithm ID : core algorithm variant
bf : brute-force
hbf : hashing-based filter
qbf : quantization-based filter and brute-force
qhbf : quantization-based filter and hashing-based filter
qpbf : quantization-based filter with pivots and brute-force

-p pivot filter mode (add -R for pivot randomization)
-n disable short-circuit
-g sequences matched in groups (for bf algorithm)
-c compact mode processing
-i interleaved bits mode (or -I with lazy evaluation)
-Q disable stats based quantization

-A ignore aligning sequences to 256-bit
-v verify results 
-q quiet output (only parameters)

```

naive brute-force algorithm benchmark (median out 9 executions) 
over Campylobacter jejuni dataset (distance limit set to 16):
```
./bench-dataset -a bf -n -r 9 campylobacter.unique.csv 16
```
benchmark using interleaved bits data representation, pivot filter 
with grouped brute-force approach over Salmonella typhi dataset (distance limit set to 89):
```
./bench-dataset -a bf -i -p -g SalmonellaTyphy_SNP.unique.txt 89
```
benchmark using brute-force algorithm and naive quantization 
(according to the lowest bit) filter with control pivot strategy 
over Campylobacter jejuni dataset (distance limit set to 8):
```
./bench-dataset -a qpbf -Q -r 9 campylobacter.unique.csv 8
```

#### Example datasets

The real datasets are freely available at [Zenodo](https://zenodo.org/record/3945108) (thanks to Alexandre P. Francisco). 

### Benchmarks on synthetic binary data - basic usage 

```
bench-binary [-a algorithmID] [-n] [-p] [-R] [-g] [-s] 
		[-r noOfRepeats] [-v] [-q] [-A] [-b bitsPerPacked] onesInPromiles m d k

algorithm ID : core algorithm variant
bf : brute-force
hbf : hashing-based filter
qbf : quantization-based filter and brute-force
qhbf : quantization-based filter and hashing-based filter
qpbf : quantization-based filter with pivots and brute-force

-p pivot filter mode (add -R for pivot randomization)
-n disable short-circuit
-g sequences matched in groups (for bf algorithm)

-A ignore aligning sequences to 256-bit
-v verify results 
-q quiet output (only parameters)

To disable binary mode apply: -b 0

```

benchmark using hash-based filter approach for generated 
random binary data with following params: density of set bits = 900 (per 1000); 
sequences length m = 512; sequences count d = 16384; Hamming distance limit k = 28
```
./bench-binary -a hbf 900 512 16384 28
```
