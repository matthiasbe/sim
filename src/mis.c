/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/mis.h"
#include <omp.h>
#include <sys/time.h>
#include <gsl/gsl_eigen.h>
#include <lapacke.h>

void estimate_error(){

}

void copy_matrix(int M, int N, double mat[M][N], double new[M][N]){
	#pragma omp parallel for
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			new[i][j] = mat[i][j];
		}
	}
}

#include <mpi.h>

void measure_accuracy(int N, int M, double q[N][M], double A[N][N], double accuracies[M], MPI_Comm comm){
	double (*transposed)[N] = (double (*)[N]) malloc(sizeof(double)*M*N);
	transpose(N, M, q, transposed);
	for (int k = 0; k < M; ++k)
	{
		double* vector = (double *) transposed[k];
		double *result = (double*) malloc(sizeof(double) * N);
		matrix_product(N, N, 1, A, (double(*)[1]) vector, (double(*)[1]) result, comm);
		double diff = scalar_product(N, vector, result)/sqrt(scalar_product(N, result, result)*scalar_product(N, vector, vector)); 
		free(result);
		accuracies[k] = diff;
	}
	
	free(transposed);
}

// Returns the scalar product of u and v
double scalar_product(int N, double u[N], double v[N]) {
	double* w = (double *) malloc(sizeof(double)*N);
	double result = 0;
	
	int i;
	#pragma omp parallel for
	for (i = 0; i < N; i++) {
		w[i] = u[i]*v[i];
	}
	
	#pragma omp parallel for reduction (+:result)
	for (i = 0; i < N; i++) {
		result = w[i] + result;
	}

	free(w);
	return result;
}

void compute_submatrix(int psize[2], int rank, int N, int M, int min[2], int max[2], MPI_Comm comm) {
	int coords[2];
	MPI_Cart_coords(comm, rank, 2, coords);
	min[0] = coords[0] * N / psize[0];
	min[1] = coords[1] * M / psize[1];
	max[0] = (coords[0] + 1) * N / psize[0] - 1;
	max[1] = (coords[1] + 1) * M / psize[1] - 1;
}

// C = AB
void matrix_product(int M, int N, int P, double A[M][N], double B[N][P], double C[M][P], MPI_Comm comm){
	
	int size;
	int psize[2], coords[2], min[2], max[2];
	
	MPI_Comm_size(comm, &size);
	MPI_Cart_coords(comm, size-1, 2, psize);
	psize[0]++;
	psize[1]++;

	MPI_Bcast(&M, 1, MPI_INT, 0, comm);
	MPI_Bcast(&N, 1, MPI_INT, 0, comm);
	MPI_Bcast(&P, 1, MPI_INT, 0, comm);
	MPI_Bcast(&(B[0][0]), N*P, MPI_DOUBLE, 0, comm);

	for(int i = 1; i<size; i++) {
		compute_submatrix(psize, i, M, P, min, max, comm);
	
		MPI_Send(&(A[min[0]][0]), (max[0] - min[0] + 1) * N, MPI_DOUBLE, i, 0, comm);
	}

	compute_submatrix(psize, 0, M, P, min, max, comm);

	double result;
	#pragma omp parallel for
	for (int k = 0; k<=max[0]; k++) {
		for (int i = 0; i<=max[1];i++) {
			result = 0;
			#pragma omp parallel for reduction (+:result)
			for (int j = 0; j<N;j++) {
				result += A[k][j] * B[j][i];
			}
			C[k][i] = result;
		}
	}

	double tmp[M*P];

	for(int i = 1; i<size; i++) {
		compute_submatrix(psize, i, M, P, min, max, comm);
		MPI_Recv(tmp, (max[0] - min[0] + 1)*(max[1]-min[1] + 1), MPI_DOUBLE, i, 2, comm, NULL);
		for(int j = 0; j <= (max[0] - min[0]); j++) {
			for (int k = 0; k <= (max[1] - min[1]); k++) {
				C[j+min[0]][k+min[1]] = tmp[j * (max[1] - min[1] + 1) + k];
			}
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
void projection(int N, int M, double A[N][N], double Z[N][M], double B[M][M], MPI_Comm comm){
	double intermediate[N][M];
	double transposed[M][N];
	matrix_product(N, N, M, A, Z, intermediate, comm);

	transpose(N, M, Z, transposed);

	matrix_product(M, N, M, transposed, intermediate, B, comm);
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
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			A[i][j] = 0;
			A[i][i] = 1 + i;
			if(j < M)
				q[i][j] = (double) rand() / RAND_MAX;
		}
	}
}

void init_q(int N, int M, double q[M][N]){
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j<N; j++) {
			q[i][j] = (double) rand() / RAND_MAX;
		}
	}
	
}

// Do the Simultaneous Iterations Methods
void mis(int N, int M, double A[N][N], double q[N][M], int iter, MPI_Comm comm) {
	
	// Temp vector
    double Z[N][M];
    double B[M][M];
    double H[M][M];
    double Zt[M][N];

    double accuracies[M];
    int return_code;
    // printf("=====A=====\n");
    // print_matrix(N, N, A);
    // printf("===========\n");
	// A^k*v serie calculation
    for(int n = 0; n < iter; n++) {
    	printf("===Accuracy===\n");
    	measure_accuracy(N, M, q, A, accuracies, comm);
        print_matrix(1, M, accuracies);
        printf("==============\n");
     //    printf("===q===\n");
     //    print_matrix(N, M, q);
     //    printf("=======\n");

		// V = A * Q
        matrix_product(N, N, M, A, q, Z, comm);

		// QR decomposition V = Z R
		transpose(N, M, Z, Zt);
        orthonormalize(M, N, Zt);
        transpose(M, N, Zt, Z);

		// B = Zt A Z
        projection(N, M, A, Z, B, comm);

        //Factorisation de Schur LAPACK
        double *tau = (double *) malloc(sizeof(double)*(M-1));
        return_code = LAPACKE_dgehrd(LAPACK_ROW_MAJOR, M, 1, M, (double *) B, M, tau);

        double (*Q)[M] = (double (*)[M]) malloc(sizeof(double)*M*M);
        double *wr = (double *) malloc(sizeof(double)*M);
        double *wi = (double *) malloc(sizeof(double)*M);
        copy_matrix(M, M, B, Q);

        return_code = LAPACKE_dorghr(LAPACK_ROW_MAJOR, M, 1, M, (double *) Q, M, tau);
        return_code = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'S', 'V', M, 1, M, (double *) B, M, (double *) wr, (double *) wi, (double *) Q, M);
        // printf("=====Q======\n");
        // print_matrix(M, M, Q);
        // printf("=============\n");

        // printf("=====B======\n");
        // print_matrix(M, M, B);
        // printf("=============\n");
        // printf("info : %d\n", return_code);
        ///////////////////////////////

        //Schur factorization B = Yt R Y
        // gsl_matrix_view gsl_B = gsl_matrix_view_array((double *)B, M, M);
        // gsl_matrix* gsl_Y = gsl_matrix_alloc(M, M);
        // gsl_vector_complex* eigenvalues = gsl_vector_complex_alloc(M);
        // gsl_eigen_nonsymm_workspace* ws = gsl_eigen_nonsymm_alloc(M);

        // gsl_eigen_nonsymm_Z(&(gsl_B.matrix), eigenvalues, gsl_Y, ws);
        // gsl_vector_view eigenvalues_real = gsl_vector_complex_real(eigenvalues);
        // gsl_vector_view eigenvalues_imag = gsl_vector_complex_imag(eigenvalues);
        for (int i = 0; i < M; ++i)
        {
        	printf("vp : %f + %fi\n", wr[i], wi[i]);
        }
        printf("\n");

		// Qk = ZY is the new approx of eigenvectors
        matrix_product(N, M, M, Z, Q, q, comm);

        // gsl_vector_complex_free(eigenvalues);
        // gsl_matrix_free(gsl_Y);
        // gsl_eigen_nonsymm_free(ws);

        // #pragma omp parallel for
        // for(int i = 0; i<M; i++) {
        //     for(int j = 0; j<N; j++) {
        //         q[i][j] = Z[i][j];
        //     }
        // }
    }
}


