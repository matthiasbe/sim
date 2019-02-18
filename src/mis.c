/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/mis.h"
#include <omp.h>
#include <sys/time.h>
#include <lapacke.h>
#include <mpi.h>
#include <float.h>

#define CSV_FILENAME "output.csv"

void swap_columns(int i, int j, int N, int M, double mat[N][M]){
	if(i != j){
		double temp;
		for (int idx = 0; idx < N; ++idx)
		{
			temp = mat[idx][i]; mat[idx][i] = mat[idx][j]; mat[idx][j] = temp;
		}
	}
}

int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

void estimate_errors(int N, int M, double Q[N][M], double Z[N][M], double eigen_real[M], double errors[M], int lock_idx){

	for (int i = lock_idx; i < M; ++i)
	{
		errors[i] = 0;
	}
	// #pragma omp parallel for simd
	for (int i = 0; i < N; ++i)
	{
		for (int j = lock_idx; j < M; ++j)
		{
			errors[j] += pow(Z[i][j] - Q[i][j]*eigen_real[j], 2);
		}
	}

}

double condition_matrix(int N, int M, double mat[N][M]){
	double max = abs(mat[0][0]);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			if(abs(mat[i][j]) > max)
				max = abs(mat[i][j]);
		}
	}

	#pragma omp parallel for simd
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			mat[i][j]/=max;
		}
	}
	return max;
}


void copy_matrix(int M, int N, double mat[M][N], double new[M][N]){
	#pragma omp parallel for simd
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			new[i][j] = mat[i][j];
		}
	}
}

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
	#pragma omp parallel for simd
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
void compute_submatrix(int psize[2], int rank, int N, int M, int min[2], int max[2], MPI_Comm comm) {

	int coords[2];
	MPI_Cart_coords(comm, rank, 2, coords);
	min[0] = coords[0] * N / psize[0];
	min[1] = coords[1] * M / psize[1];
	max[0] = (coords[0] + 1) * N / psize[0] - 1;
	max[1] = (coords[1] + 1) * M / psize[1] - 1;
}

void matrix_product_wlock(int M, int N, int P, double A[M][N], double B[N][P], double C[M][P], int lock_idx){
	double result;
	#pragma omp parallel for
	for (int k = 0; k<M; k++) {
		for (int i = 0; i < lock_idx; ++i)
		{
			C[k][i] = B[k][i];
		}
		for (int i = lock_idx; i<P;i++) {
			result = 0;
			#pragma omp parallel for reduction (+:result)
			for (int j = 0; j<N;j++) {
				result += A[k][j] * B[j][i];
			}
			C[k][i] = result;
		}
	}
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
	#pragma omp parallel for simd
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			T[i][j] = A[j][i];
		}
	}
}

void crop_columns(int N, int M, int krylov_size, int lock_idx, double mat[N][M], double cropped[N][krylov_size]){
	for (int i = 0; i < N; ++i)
	{
		for (int j = lock_idx; j < M; ++j)
		{
			cropped[i][j-lock_idx] = mat[i][j];
		}
	}
}

void merge_in(int N, int M, int krylov_size, double dst[N][M], double src[N][krylov_size]){
	for (int i = 0; i < N; ++i)
	{
		for (int j = M - krylov_size; j < M; ++j)
		{
			dst[i][j] = src[i][j - M + krylov_size];
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

/* Computes B=q(^H)Aq
*/
void projection_wlock(int N, int krylov_size, double A[N][N], double Zcropped[N][krylov_size], double B[krylov_size][krylov_size], int lock_idx){
	double intermediate[N][krylov_size];
	double transposed[krylov_size][N];

	transpose(N, krylov_size, Zcropped, transposed);

	matrix_product_wlock(N, N, krylov_size, A, Zcropped, intermediate, 0);
	matrix_product_wlock(krylov_size, N, krylov_size, transposed, intermediate, B, 0);
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
			printf("[%e]",A[i][j]);
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
			A[i][j] = (double) rand() / RAND_MAX;
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
void mis(int N, int M, double A[N][N], double q[N][M], int iter, int precision, int n_eigen, MPI_Comm comm) {
	// double norm_factor = condition_matrix(N, N, A);
	// Temp vector
    double Z[N][M];
    double Zt[M][N];
    double eigen_real[M];

    double accuracies[M];
    int return_code;
    int lock_idx = 0;
    int krylov_size = M;

	struct timeval start, curr;
	gettimeofday(&start, NULL);
	
	// This limit iterations to 150 in case only a precision is given
	if (iter == 0) iter = 10000;

	FILE* fp = fopen(CSV_FILENAME, "w+");
	fprintf(fp, "iteration, eigindex, value, type\n");
	fprintf(fp, "0,N/A,%d,t\n", 0);
	
	// V = A * Q
    matrix_product(N, N, M, A, q, Z, comm);

    for(int n = 0; n < iter; n++) {
     	// measure_accuracy(N, M, q, A, accuracies, comm);

		// QR decomposition V = Z R
		transpose(N, M, Z, Zt);
        orthonormalize(M, N, Zt);
        transpose(M, N, Zt, Z);
        
        // printf("=====Z ORTHONORMALIZED=====\n");
        // print_matrix(N, M, Z);
        // printf("===========================\n");


        double B[krylov_size][krylov_size];
   		double H[krylov_size][krylov_size];
   		double Zcropped[N][krylov_size];
		// B = Zt A Z
        // projection(N, M, A, Z, B, comm);
        crop_columns(N, M, krylov_size, lock_idx, Z, Zcropped);
        // printf("=====Z CROPPED=====\n");
        // print_matrix(N, krylov_size, Zcropped);
        // printf("===========================\n");
        projection_wlock(N, krylov_size, A, Zcropped, B, lock_idx);
        // printf("=====B=====\n");
        // print_matrix(krylov_size, krylov_size, B);
        // printf("===========================\n");

        //Factorisation de Schur LAPACK B = Yt R Y
        double *tau = (double *) malloc(sizeof(double)*(krylov_size-1));
        return_code = LAPACKE_dgehrd(LAPACK_ROW_MAJOR, krylov_size, 1, krylov_size, (double *) B, krylov_size, tau);

        double (*Q)[krylov_size] = (double (*)[krylov_size]) malloc(sizeof(double)*krylov_size*krylov_size);
        double *wi = (double *) malloc(sizeof(double)*krylov_size);
        copy_matrix(krylov_size, krylov_size, B, Q);

        return_code = LAPACKE_dorghr(LAPACK_ROW_MAJOR, krylov_size, 1, krylov_size, (double *) Q, krylov_size, tau);
        return_code = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'S', 'V', krylov_size, 1, krylov_size, (double *) B, krylov_size, ((double *) eigen_real) + lock_idx, (double *) wi, (double *) Q, krylov_size);

        for (int i = 0; i < krylov_size; ++i)
        {
			fprintf(fp,"%d,%d,%f,vr\n", n+1, i + lock_idx, eigen_real[i + lock_idx]);
			fprintf(fp,"%d,%d,%f,vi\n", n+1, i + lock_idx, wi[i]);
        }

		// Qk = ZY is the new approx of eigenvectors
		double qcropped[N][krylov_size];
        matrix_product(N, krylov_size, krylov_size, Zcropped, Q, qcropped, comm);
        // printf("=====q CROPPED=====\n");
        // print_matrix(N, krylov_size, qcropped);
        // printf("===========================\n");
        merge_in(N, M, krylov_size, q, qcropped);
        // printf("=====q MERGED=====\n");
        // print_matrix(N, M, q);
        // printf("===========================\n");
        free(tau);
        free(wi);
        free(Q);

		gettimeofday(&curr, NULL);
		double duration = (double) (curr.tv_usec - start.tv_usec) / 1000000 +
		         (double) (curr.tv_sec - start.tv_sec);
		fprintf(fp, "%d,N/A,%f,t\n", n+1, duration);

		// V = A * Q
        matrix_product_wlock(N, N, M, A, q, Z, lock_idx);
        // printf("=====NEW Z=====\n");
        // print_matrix(N, M, Z);
        // printf("===========================\n");

		estimate_errors(N, M, q, Z, eigen_real, accuracies, lock_idx);
		// print_matrix(1, M, eigen_real);
		double accuracies_sorted[M];
		for (int i = 0; i < M; ++i)
		{
			accuracies_sorted[i] = accuracies[i];
		}
    	qsort(accuracies_sorted, M, sizeof(double), compare_doubles);
    	// print_matrix(1, M, accuracies);
		for (int i = 0; i < n_eigen; i++) {
			fprintf(fp, "%d,%d,%e,A\n", n, i, accuracies_sorted[i]);
		}
		if (precision > 0){
			if(accuracies_sorted[n_eigen - 1] < pow(10,-precision)) {
				printf("**** accuracy %d reached with ****\n", precision);
				printf("minimum eigenvector precision : %f\n", accuracies_sorted[n_eigen - 1]);
				printf("Number of iteration : %d\n", n);
				break;
			}
			else{
				for (int i = lock_idx; i < M; ++i)
				{
					if (accuracies[i] < pow(10,-precision))
					{
						printf("Iteration : %d => SWAPPING %d and %d because accuracy is %e\n", n, lock_idx, i, accuracies[i]);
						swap_columns(lock_idx, i, N, M, Z);
						swap_columns(lock_idx, i, N, M, q);
						double temp = accuracies[i]; accuracies[i] = accuracies[lock_idx]; accuracies[lock_idx] = temp;
						temp = eigen_real[i]; eigen_real[i] = eigen_real[lock_idx]; eigen_real[lock_idx] = temp;
						++lock_idx;
						--krylov_size;
					}
				}
			}
		}

		if (iter < 20 || n % (iter / 20) == 0)
			printf("%2d%% (%d iterations) in %.3f seconds\n", (n*100)/iter, n, duration);
    }
}


