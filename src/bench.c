#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include "mis.h"
#include "read_file.h"
#include <math.h>

#include <mpi.h>

#define MPI_TERM_CODE -1

struct arguments {
	// Size of the A matrix
	int N;
	// Number of eigenvectors to guess
	int M;
	// Size of Krylov subspace
	int e;
	// Number of eigenvectors to compute
	int iter;
	// A matrix passed as a file
	char* matrix_filename;
};

void parse_args(struct arguments* args, int argc, char* argv[]) {
	char c;

	args->N = 0;
	args->M = 0;
	args->e = 0;
	args->iter = 0;
	args->matrix_filename = NULL;

	while ((c = getopt (argc, argv, "e:n:m:i:f:h")) != -1) {
		switch (c)
		{
			case 'n':
				args->N = atoi(optarg);
				break;
			case 'm':
				args->M = atoi(optarg);
				break;
			case 'e':
				args->e = atoi(optarg);
				break;
			case 'i':
				args->iter = atoi(optarg);
				break;
			case 'f':
				args->matrix_filename = optarg;
				break;
			case 'h':

				printf("Usage : bench -n <matrix-size> -m <Krylov_size> -i <iter-nb> -e <eigvec-nb>\n");
				printf("Or  	bench -i <iter-nb> -f <matrix-filename> -e <eigvec-nb>\n");
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

void run_worker(int rank, int size, MPI_Comm comm) {
	// C matrix coordinate to consider for computing C = A*B
	// min[0] first line index
	// max[0] last line index
	// min[1] first column index
	// max[1] last column index
	
	//               (  C global matrix   )
	//               (                    )
	//       min[0] -( - - +------+       )
	//               (     |local |       )
	//               (     |C     |       )
	//       max[0] -( - - +------+       )
	//               (                    )
	//               (     |      |       )
	//					min[1]   max[1]
	int pcoord[2], min[2], max[2], psize[2];

	// Total matrices size : A = MxN  -  B = NxP  -   C = MxP
	// Local matrices size :
	// 		A = (max[0] - min[0]) x N
	// 		B = N x P
	// 		C = (max[0] - min[0]) x (max[1] - min[1])
	// Considered matrices for local computating
	// 		A = (max[0] - min[0]) x N
	// 		B = N x (max[1] - min[1])
	// 		C = (max[0] - min[0]) x (max[1] - min[1])
	
	// M can be -1 meaning main process is over
	int M, N, P;

	// Get pcoord = coord of current proc in process grid
	MPI_Cart_coords(comm, rank, 2, pcoord);

	// Get psize = dimension of process grid
	MPI_Cart_coords(comm, size-1, 2, psize);
	psize[0]++;
	psize[1]++;

	MPI_Bcast(&M, 1, MPI_INT, 0, comm);

	while(M != MPI_TERM_CODE) {

		MPI_Bcast(&N, 1, MPI_INT, 0, comm);
		MPI_Bcast(&P, 1, MPI_INT, 0, comm);

		// Get min and max
		compute_submatrix(psize, rank, M, P, min, max, comm);

		double *A = malloc((max[0] - min[0] + 1)*N*sizeof(double));
		double *B = malloc(N*P*sizeof(double));
		double *C = malloc((max[0] - min[0] + 1)*(max[1] - min[1] + 1)*sizeof(double));

		// Receive matrices
		MPI_Bcast(B, N*P, MPI_DOUBLE, 0, comm);
		MPI_Recv(A, (max[0] - min[0] + 1) * N, MPI_DOUBLE, 0, 0, comm, NULL);

		// Compute C = A*B[:][min[1] , max[1]]
		double result;
		#pragma omp parallel for
		for (int k = 0; k<=max[0] - min[0]; k++) {
			for (int i = min[1]; i<=max[1];i++) {
				result = 0;
				#pragma omp parallel for reduction (+:result)
				for (int j = 0; j<N;j++) {
					result += A[k*N+j] * B[j*P+i];
				}
				C[k*(max[1] - min[1] + 1)+i - min[1]] = result;
			}
		}

		MPI_Send(C, (max[0] - min[0] + 1)*(max[1]-min[1] + 1), MPI_DOUBLE, 0, 2, comm);

		free(A);
		free(B);
		free(C);

		MPI_Bcast(&M, 1, MPI_INT, 0, comm);
	}
}

void run_main_process(struct arguments args, MPI_Comm comm) {
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
		A = (double (*)[args.N]) malloc(sizeof(double)*args.N*args.N);
		q = (double (*)[args.M]) malloc(sizeof(double)*args.M*args.N);
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

	int term_code = MPI_TERM_CODE;
	MPI_Bcast(&term_code, 1, MPI_INT, 0, comm);
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


	if (rank == 0) {
		
		struct arguments args;
		parse_args(&args, argc, argv);


		run_main_process(args, comm);
	}
	else {
		run_worker(rank, size, comm);
	}

	MPI_Finalize();
}
