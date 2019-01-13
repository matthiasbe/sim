#! /bin/sh

gcc -fopenmp -g -Wall mis.c test.c read_file.c -o test -lm
gcc -fopenmp -g -Wall mis.c bench.c -o bench -lm
