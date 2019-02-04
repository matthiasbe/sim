#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include "mis.h"
#include "read_file.h"

struct arguments {
	// Size of the A matrix
	int N;
	// Number of eigenvectors to guess
	int M;
	// Number of iterations
	int iter;
	// A matrix passed as a file
	char* matrix_filename;
};

void parse_args(struct arguments* args, int argc, char* argv[]) {
	char c;

	args->N = 0;
	args->M = 0;
	args->iter = 0;
	args->matrix_filename = NULL;

	while ((c = getopt (argc, argv, "n:m:i:f:h")) != -1) {
		switch (c)
		{
			case 'n':
				args->N = atoi(optarg);
				break;
			case 'm':
				args->M = atoi(optarg);
				break;
			case 'i':
				args->iter = atoi(optarg);
				break;
			case 'f':
				args->matrix_filename = optarg;
				break;
			case 'h':
				printf("Usage : bench -n <matrix-size> -m <eigvec-nb> -i <iter-nb>\n");
				printf("Or  	bench -i <iter-nb> -f <matrix-filename>\n");
				exit(0);
			case '?':
				if (optopt == 'n' || optopt == 'm' || optopt == 'i')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr,
							"Unknown option character `\\x%x'.\n",
							optopt);
				printf("Use -h to print help\n");
				exit(-1);
			default:
				abort ();
		}
	}

	if (args->iter == 0 || args->M == 0 ||
			((args->matrix_filename == NULL) == (args->N == 0))) {
		fprintf(stderr, "bad arguments\n");
		printf("Use -h option for help\n");
		exit(-1);
	}
}

int main(int argc, char* argv[]) {

	struct arguments args;

	parse_args(&args, argc, argv);

	double (*A)[args.N] = (double (*)[args.N]) malloc(sizeof(double)*args.N*args.N);
	double (*q)[args.N] = (double (*)[args.N]) malloc(sizeof(double)*args.M*args.N);

	init(args.N, args.M, A, q);

	struct timeval start;
	gettimeofday(&start, NULL);

	mis(args.N, args.M, A, q, args.iter);

	struct timeval end;
	gettimeofday(&end, NULL);
	double duration = (double) (end.tv_usec - start.tv_usec) / 1000000 +
		         (double) (end.tv_sec - start.tv_sec);
	printf("duration (s) : %f\n", duration);

	free(A);
	free(q);
}
