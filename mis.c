/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Size of the A matrix
#define N 5000
// Number of eigenvectors to guess
#define M 50
// Number of iterations
#define ITER 50

// Returns the scalar product of u and v
double scalar_product(double u[N], double v[N]) {

	double* w = malloc(sizeof(double)*N);
	double result = 0;
	
	int i;
	#pragma omp parallel for num_threads(4)
	for (i = 0; i < N; i++) {
		w[i] = u[i]*v[i];
	}
	
	#pragma omp parallel for reduction (+:result) num_threads(4)
	for (i = 0; i < N; i++) {
		result = w[i] + result;
	}

	free(w);
	return result;
}

/* Make the system of column vecotrs of the matrix A orthogonal and orthonormal
 * using Gram-Schmidt algorithm */
void orthonormalize(double A[M][N]) {
	
	double temp[N];
	double norm, temp_norm;

	for(int i = 0; i<M; i++) {
		for (int j = 0; j<N; j++) {
			temp[j] = A[i][j];
		}

		for (int k = 0; k<i; k++) {
			norm = scalar_product(A[i], A[k]);
			norm /= scalar_product(A[k], A[k]);


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
void print_matrix(double A[N][N]) {
	for (int i = 0; i<N ; i++) {
		for (int j = 0; j<N; j++) {
			printf("[%f]",A[i][j]);
		}
		printf("\n");
	}
}

// Initialize A and q with random values
void init(double *A[N], double *q[M]) {
	for (int i = 0; i<N; i++) {
		A[i] = malloc(sizeof(double) * N);
		if (i < M) {
			q[i] = malloc(sizeof(double) * N);
		}
	}
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
void mis(double *A[N], double *q[M]) {
	// Temp vector
        double v[M][N];

	// A^k*v serie calculation
        for(int n = 0; n < ITER; n++) {

		// v = A * Q
		for (int k = 0; k<M; k++) {
			for (int i = 0; i<N;i++) {
				v[k][i] = 0.0;
				for (int j = 0; j<3;j++) {
					v[k][i] += q[k][j]*A[i][j];
				}
				v[k][i] = v[k][i];
			}
		}
		
		orthonormalize(v);

		// q = v
                for(int i = 0; i<M; i++) {
			for(int j = 0; j<N; j++) {
				q[i][j] = v[i][j];
			}
                }
		
        }
}

int main() {

	double *A[N];
        double *q[M];

	init(A, q);

	mis(A, q);
}


