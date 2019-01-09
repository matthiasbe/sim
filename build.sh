#! /bin/sh

gcc -fopenmp -g -Wall mis.c -o mis -lm -llapacke -llapack -I/usr/include/lapacke -L/usr/lib64/
