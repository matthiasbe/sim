#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include "mis.h"
#include "read_file.h"
#include <math.h>

#include <mpi.h>

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

	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm comm;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int psize[2];
	int cart_period[2] = { 0, 0 };

	// number of processes in the x dimension
	psize[0] = sqrt(size);
	if ( psize[0]<1 ) psize[0]=1; 
	// number of processes in the y dimension
	psize[1] = size / psize[0];
	if (psize[0]*psize[1] != size) {
        fprintf(stderr, "Error: invalid number of processes\n");
        abort();
    }

	MPI_Cart_create(MPI_COMM_WORLD, 2, psize, cart_period, 1, &comm);
	MPI_Comm_rank(comm, &rank);

	printf("size : %i (%i x %i)\n", size, psize[0], psize[1]);

	printf("Teeeest : %d\n", rank);

	if (rank == 0) {
		struct arguments args;

		parse_args(&args, argc, argv);


		double *mat;
		double (*A)[args.N];
		double (*q)[args.N];

		if (args.matrix_filename != NULL) {
			int size[2];
			read_mtx(args.matrix_filename, size, &mat);
			if (size[0] != size[1]) {
				fprintf(stderr, "Cannot handle non-square matrices.\n");
				exit(-1);
			}

			args.N = size[0];
			A = (double (*) []) mat;
			q = (double (*) []) malloc(sizeof(double)*args.M*args.N);
			init_q(args.N, args.M, q);

		} else {
			double (*A)[args.N] = (double (*)[args.N]) malloc(sizeof(double)*args.N*args.N);
			double (*q)[args.N] = (double (*)[args.N]) malloc(sizeof(double)*args.M*args.N);
			init(args.N, args.M, A, q);
		}

		struct timeval start;
		gettimeofday(&start, NULL);

		mis(args.N, args.M, A, q, args.iter, comm);

		struct timeval end;
		gettimeofday(&end, NULL);
		double duration = (double) (end.tv_usec - start.tv_usec) / 1000000 +
					 (double) (end.tv_sec - start.tv_sec);
		printf("duration (s) : %f\n", duration);

		free(A);
		free(q);

		MPI_Abort(MPI_COMM_WORLD, 0);
	}
	else {
		int pcoord[2], min[2], max[2], psize[2];
		int M, N, P;
		printf("worker : %i\n", rank);
		MPI_Cart_coords(comm, rank, 2, pcoord);
		MPI_Cart_coords(comm, size-1, 2, psize);
		psize[0]++;
		psize[1]++;

		while(1) {
			MPI_Recv(&M, 1, MPI_INT, 0, 0, comm, NULL);
			MPI_Recv(&N, 1, MPI_INT, 0, 0, comm, NULL);
			MPI_Recv(&P, 1, MPI_INT, 0, 0, comm, NULL);

			compute_submatrix(psize, rank, M, P, min, max, comm);

			double *A = malloc((max[0] - min[0])*N*sizeof(double));
			double *B = malloc(N*P*sizeof(double));
			double *C = malloc((max[0] - min[0])*(max[1] - min[1])*sizeof(double));
			printf("[%i] received M=%i N=%i P=%i\n", rank, M, N, P);

			MPI_Recv(A, (max[0] - min[0]) * N, MPI_DOUBLE, 0, 0, comm, NULL);
			printf("[%i]received A\n", rank);
			MPI_Recv(B, N*P, MPI_DOUBLE, 0, 1, comm, NULL);
			printf("[%i]received B\n", rank);

//
//			double result;
//			#pragma omp parallel for
//			for (int k = 0; k<max[0] - min[0]; k++) {
//				for (int i = min[1]; i<max[1];i++) {
//					result = 0;
//					#pragma omp parallel for reduction (+:result)
//					for (int j = 0; j<N;j++) {
//						result += A[k][j] * B[j][i];
//					}
//					C[k][i] = result;
//				}
//			}

			MPI_Send(C, (max[0] - min[0])*(max[1]-min[1]), MPI_DOUBLE, 0, 2, comm);
		}

	}

	MPI_Finalize();
}
