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
double scalar_product(int N, double u[N], double v[N]) {
	double* w = (double *) malloc(sizeof(double)*N);
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
	//#pragma omp parallel for
	for (int k = 0; k<M; k++) {
		for (int i = 0; i<P;i++) {
			result = 0.0;
			//#pragma omp parallel for reduction (+:result)
			for (int j = 0; j<N;j++) {
				result += A[k][j] * B[j][i];
			}
			C[k][i] = result;
		}
	}
}

void transpose(int M, int N, double A[M][N], double T[N][M]){
	//#pragma omp parallel for
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

/* Make the system of line vectors of the matrix A orthogonal and orthonormal
 * using Gram-Schmidt algorithm */
void orthonormalize(int N, int M, double A[N][M]) {
	double temp_norm;
	for(int i = 1; i<N; i++) {
		for (int k = 0; k<i; k++) {
			double norm = scalar_product(M, A[i], A[k]) / scalar_product(M, A[k], A[k]);
			for (int j = 0; j<M; j++) {
				A[i][j] -= norm * A[k][j];
			}
		}
	}

	for(int i = 0; i<N; i++) {
		temp_norm = 0;
		for (int k = 0; k<M; k++) {
			temp_norm += A[i][k] * A[i][k];
		}

		temp_norm = sqrt(temp_norm);

		for (int j = 0; j<M; j++) {
			A[i][j] /= temp_norm;
		}
	}
}

// Print the A matrix for debug purpose
void print_matrix(int N, int M, double A[N][M]) {
	for (int i = 0; i<N ; i++) {
		for (int j = 0; j<M; j++) {
			printf("[%f]",A[i][j]);
		}
		printf("\n");
	}
}

// Initialize A and q with random values
void init(int N, int M, double A[N][N], double q[N][M]) {
	for (int i = 0; i<N; i++) {
		double random_dbl = (double) rand() / RAND_MAX;
		for (int j = 0; j<i; j++) {
			A[i][j] = random_dbl;
			A[j][i] = random_dbl;
			random_dbl = (double) rand() / RAND_MAX;
		}
		A[i][i] = random_dbl;

		for (int j = 0; j<M; j++) {
			random_dbl = (double) rand() / RAND_MAX;
			q[i][j] = random_dbl;
		}
	}
}

// Do the Simultaneous Iterations Methods
void mis(int N, int M, double A[N][N], double q[N][M], int iter) {
	// Temp vector
    double Z[N][M];
    double B[M][M];
    double H[M][M];
    double Zt[M][N];

	// A^k*v serie calculation
    for(int n = 0; n < iter; n++) {
    	// printf("Iteration number %d\n", n);
    	// printf("A =\n");
    	// print_matrix(N, N, A);
    	// printf("\n");

    	// printf("q =\n");
    	// print_matrix(N, M, q);
    	// printf("\n");

		// V = A * Q
        matrix_product(N, N, M, A, q, Z);

     //    printf("V =\n");
    	// print_matrix(N, M, Z);
    	// printf("\n");
	
		// QR decomposition V = Z R
		transpose(N, M, Z, Zt);
        orthonormalize(M, N, Zt);
        transpose(M, N, Zt, Z);

     //    printf("Z =\n");
    	// print_matrix(N, M, Z);
    	// printf("\n");

		// B = Zt A Z
        projection(N, M, A, Z, B);

     //    printf("B =\n");
    	// print_matrix(M, M, B);
    	// printf("\n");

        //Schur factorization B = Yt R Y
        gsl_matrix_view gsl_B = gsl_matrix_view_array((double *)B, M, M);
        gsl_matrix* gsl_Y = gsl_matrix_alloc(M, M);
        gsl_vector_complex* eigenvalues = gsl_vector_complex_alloc(M);
        gsl_eigen_nonsymm_workspace* ws = gsl_eigen_nonsymm_alloc(M);

        gsl_eigen_nonsymm_Z(&(gsl_B.matrix), eigenvalues, gsl_Y, ws);

        // printf("eigenvalues =\n");
    	// print_matrix(1, M, eigenvalues->data);
    	// printf("\n");

    	// printf("Y =\n");
    	// print_matrix(M, M, gsl_Y->data);
    	// printf("\n");


		// Qk = ZY is the new approx of eigenvectors
        matrix_product(N, M, M, Z, (double (*)[M]) gsl_Y->data, q);

        gsl_vector_complex_free(eigenvalues);
        gsl_matrix_free(gsl_Y);
        gsl_eigen_nonsymm_free(ws);

        #pragma omp parallel for
        for(int i = 0; i<M; i++) {
            for(int j = 0; j<N; j++) {
                q[i][j] = Z[i][j];
            }
        }
    }
}


