#! /bin/sh

gcc -llapacke -fopenmp -g -Wall mis.c -o mis -lm -llapack -lgfortran
