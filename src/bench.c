#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mis.h"

int main(int argc, char* argv[]) {

	if (argc != 4) {
		printf("Usage : ./bench N M <nb-iterations>\n");
		printf("N : size of the matrix\n");
		printf("M : number of eigenvectors to guess\n");
		exit(-1);
	}

	// Size of the A matrix
	int N = atoi(argv[1]);
	// Number of eigenvectors to guess
	int M = atoi(argv[2]);
	// Number of iterations
	int iter = atoi(argv[3]);

	double (*A)[N] = (double (*)[N]) malloc(sizeof(double)*N*N);
	double (*q)[N] = (double (*)[N]) malloc(sizeof(double)*M*N);

	init(N, M, A, q);

	struct timeval start;
	gettimeofday(&start, NULL);

	mis(N, M, A, q, iter);

	struct timeval end;
	gettimeofday(&end, NULL);
	double duration = (double) (end.tv_usec - start.tv_usec) / 1000000 +
		         (double) (end.tv_sec - start.tv_sec);
	printf("duration (s) : %f\n", duration);

	free(A);
	free(q);
}
