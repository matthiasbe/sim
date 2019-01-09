/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
// #include <lapacke.h>

// Returns the scalar product of u and v
double scalar_product(int N, int M, double u[N], double v[N]) {

	double* w = malloc(sizeof(double)*N);
	double result = 0;
	
	int i;
	//#pragma omp parallel for
	for (i = 0; i < N; i++) {
		w[i] = u[i]*v[i];
	}
	
	//#pragma omp parallel for reduction (+:result)
	for (i = 0; i < N; i++) {
		result = w[i] + result;
	}

	free(w);
	return result;
}

/* Make the system of column vecotrs of the matrix A orthogonal and orthonormal
 * using Gram-Schmidt algorithm */
void orthonormalize(int N, int M, double A[M][N]) {
	
	double temp[N];
	double norm, temp_norm;

	for(int i = 0; i<M; i++) {
		for (int j = 0; j<N; j++) {
			temp[j] = A[i][j];
		}

		for (int k = 0; k<i; k++) {
			norm = scalar_product(N, M, A[i], A[k]);
			norm /= scalar_product(N, M, A[k], A[k]);


			for (int j = 0; j<N; j++) {
				temp[j] -= norm * A[k][j];
			}
		}
		for (int j = 0; j<N; j++) {
			A[i][j] = temp[j];
		}
	}

	for(int i = 0; i<M; i++) {

		temp_norm = 0;
		for ( int k = 0; k<N; k++) {
			temp_norm += A[i][k] * A[i][k];
		}

		temp_norm = sqrt(temp_norm);

		for (int j = 0; j<N; j++) {
			A[i][j] /= temp_norm;
		}
	}


}

// Print the A matrix for debug purpose
void print_matrix(int N, double A[N][N]) {
	for (int i = 0; i<N ; i++) {
		for (int j = 0; j<N; j++) {
			printf("[%f]",A[i][j]);
		}
		printf("\n");
	}
}

// Initialize A and q with random values
void init(int N, int M, double (*A)[N], double (*q)[N]) {
	for (int i = 0; i<N; i++) {
		double random_dbl = (double) rand();
		for (int j = 1; j<i; j++) {
			A[i][j] = random_dbl;
			A[j][i] = random_dbl;
			random_dbl = (double) rand();
		}
		A[i][i] = random_dbl;

		for (int j = 0; j<M; j++) {
			random_dbl = (double) rand();
			q[j][i] = random_dbl;
		}
	}
}

// Do the Simultaneous Iterations Methods
void mis(int N, int M, double (*A)[N], double (*q)[M], int iter) {
	// Temp vector
        double v[M][N];
	double result;

	// A^k*v serie calculation
        for(int n = 0; n < iter; n++) {

		// v = A * Q
		#pragma omp parallel for
		for (int k = 0; k<M; k++) {
			//#pragma omp parallel for
			for (int i = 0; i<N;i++) {
				result = 0.0;
				//#pragma omp parallel for reduction (+:result)
				for (int j = 0; j<N;j++) {
					result += q[k][j]*A[i][j];
				}
				v[k][i] = result;
			}
		}
		
		orthonormalize(N, M, v);

		// q = v
		#pragma omp parallel for
                for(int i = 0; i<M; i++) {
			for(int j = 0; j<N; j++) {
				q[i][j] = v[i][j];
			}
                }
		
        }
}


// void check_egeinvectors(int N, int M, double *A[N], double *q[N]) {
// 	double D[N];
// 	double E[N-1];
// 	double TAU[N-1];
// 	double WORK[1];
// 	int LWORK = 1;
// 	int INFO;
// 	dsytrd_("U", &N, (double *)A, &N, D, E, (double *)TAU, WORK, &LWORK, &INFO);

// 	double precision = 1;
// 	int NSPLIT;
// 	double *eigenvalues = malloc(N*sizeof(double));
//        	int IBLOCK[N], ISPLIT[N];
// 	int P = N-M;
// 	dstebz_("I", "E", &N, 0, 0, &P, &N, &precision, D, E, &M, &NSPLIT, eigenvalues, IBLOCK, ISPLIT, WORK, &LWORK, &INFO);

// 	if (INFO != 0) {
// 		printf("Error finding EV : %d\n", INFO);
// 	}
	
// 	for ( int i; i<N ; i++) {
// 		printf("[%f]", eigenvalues[i]);
// 	}
// 	printf("\n");
// }

int main(int argc, char* argv[]) {

	if (argc != 4) {
		printf("Usage : ./mis N M <nb-iterations>\n");
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

	mis(N, M, A, q, iter);

	free(A);
	free(q);
	// check_egeinvectors(N, M, A, q);
}


