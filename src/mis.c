/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include <gsl/gsl_eigen.h>

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
// C = AB
void matrix_product(int M, int N, int P, double A[M][N], double B[N][P], double C[M][P]){
	double result;
	#pragma omp parallel for
	for (int k = 0; k<M; k++) {
		for (int i = 0; i<P;i++) {
			result = 0.0;
			#pragma omp parallel for reduction (+:result)
			for (int j = 0; j<N;j++) {
				result += A[k][j] * B[i][j];
			}
			C[k][i] = result;
		}
	}
}

void transpose(int M, int N, double A[M][N], double T[N][M]){
	#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			T[i][j] = A[j][i];
		}
	}
}

/* Computes B=q(^H)Aq
*/
void projection(int N, int M, double A[N][N], double Z[N][M], double B[M][M]){
	double intermediate[N][M];
	double transposed[M][N];
	matrix_product(N, N, M, A, Z, intermediate);

	transpose(N, M, Z, transposed);

	matrix_product(M, N, M, transposed, intermediate, B);
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
void init(int N, int M, double A[N][N], double q[N][M]) {
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
void mis(int N, int M, double A[N][N], double q[N][M], int iter) {
	// Temp vector
    double Z[N][M];
    double B[M][M];
    double H[M][M];

	// A^k*v serie calculation
    for(int n = 0; n < iter; n++) {
		// v = A * Q
        matrix_product(N, N, M, A, q, Z);
	
	orthonormalize(N, M, Z);

	projection(N, M, A, Z, B);

	//Schur factorization
	gsl_matrix_view gsl_B = gsl_matrix_view_array((double *)B, M, M);
	gsl_matrix* gsl_Y = gsl_matrix_alloc(M, M);
	gsl_vector_complex* eigenvalues = gsl_vector_complex_alloc(M);
	gsl_eigen_nonsymm_workspace* ws = gsl_eigen_nonsymm_alloc(M);

	gsl_eigen_nonsymm_Z(&(gsl_B.matrix), eigenvalues, gsl_Y, ws);


	matrix_product(N, M, M, Z, (double (*)[M]) gsl_Y->data, q);

	gsl_vector_complex_free(eigenvalues);
	gsl_matrix_free(gsl_Y);
	gsl_eigen_nonsymm_free(ws);

	// q = v
	#pragma omp parallel for
	for(int i = 0; i<M; i++) {
		for(int j = 0; j<N; j++) {
			q[i][j] = Z[i][j];
		}
	}
    }
}

// int main(int argc, char* argv[]) {

// 	if (argc != 4) {
// 		printf("Usage : ./mis N M <nb-iterations>\n");
// 		printf("N : size of the matrix\n");
// 		printf("M : number of eigenvectors to guess\n");
// 		exit(-1);
// 	}

// 	// Size of the A matrix
// 	int N = atoi(argv[1]);
// 	// Number of eigenvectors to guess
// 	int M = atoi(argv[2]);
// 	// Number of iterations
// 	int iter = atoi(argv[3]);

// 	double (*A)[N] = (double (*)[N]) malloc(sizeof(double)*N*N);
//     double (*q)[N] = (double (*)[N]) malloc(sizeof(double)*M*N);

// 	init(N, M, A, q);

// 	struct timeval start;
// 	gettimeofday(&start, NULL);

// 	mis(N, M, A, q, iter);

// 	struct timeval end;
// 	gettimeofday(&end, NULL);
// 	double duration = (double) (end.tv_usec - start.tv_usec) / 1000000 +
// 		         (double) (end.tv_sec - start.tv_sec);
// 	printf("duration (s) : %f\n", duration);

// 	free(A);
// 	free(q);
// }


