#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../src/mis.h"
#include "../src/read_file.h"
#include <limits.h>

void init_q(int N, int M, double q[M][N]){
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j<N; j++) {
			q[i][j] = (double) rand();
		}
	}
	
}

void normalize(int N, double vect[N]){
	double norm = sqrt(scalar_product(N, vect, vect));
	for (int i = 0; i < N; ++i)
	{
		vect[i] /= norm;
	}
}

void test3x3() {
	int size[2];
	double *mat;
	read_matrix("../test/matrices/simple3x3", size, &mat);
	double (*A)[size[1]] = (double (*) []) mat;

	//double expected_eval[3] = {232097.75470006, -231154.01191683, -778.12519323};

	//double expected_evec[3][3] = {
    //    {0.576916252, 0.577988900, 0.0034060215},
    //    {0.577580793, -0.577016915, 0.0999994065},
    //    {0.577553517, -0.577044462, 0.000517927165}
	//};
    
    double expected_evec[3][3] = {
        { 0.34691808,  0.93639405,  0.05304737 },
        { 0.71085034, -0.22561874, -0.66617414 },
        { 0.61183302, -0.26881659,  0.74390725 }
    };

	double computed_evec[3][3];
	init_q(3, 3, computed_evec);

		mis(3, 3, (double (*)[]) A, computed_evec, 250);

		for (int i = 0; i<3; i++) {
			for ( int j = 0; j<3; j++) {
				printf("%f, ", computed_evec[i][j] / expected_evec[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	free(A);
}

void test_eigenvector(int N, int M){
	printf("Test eigenvector computation of size %dx%d\n", N, N);
	int size[2];
	// read_matrix(matrix_filename, size, &mat);
	double (*A)[N] = (double (*)[N]) malloc(sizeof(double)*N*N);
	double (*q)[M] = (double (*)[M]) malloc(sizeof(double)*M*N);
	init(N, M, A, q);
	// print_matrix(N, N, A);

	mis(N, M, A, q, 10000);

	double (*transposed)[N] = (double (*)[N]) malloc(sizeof(double)*M*N);
	transpose(N, M, q, transposed);
	free(q);
	for (int k = 0; k < M; ++k)
	{
		double* vector = (double *) transposed[k];
		double *result = (double*) malloc(sizeof(double) * N);
		matrix_product(N, N, 1, A, (double(*)[1]) vector, (double(*)[1]) result);
		printf("Scalar prod 1 : %lf\n", scalar_product(N, vector, result));
		printf("Scalar prod 2 : %lf\n", scalar_product(N, result, result));
		printf("Scalar prod 3 : %lf\n", scalar_product(N, vector, vector));
		double diff = scalar_product(N, vector, result)/(scalar_product(N, result, result)*scalar_product(N, vector, vector)); 
		free(result);
		printf("Got difference of %lf for eigenvector %d\n", diff, k);
	}
	
	free(A);
	free(transposed);
	
}

void test_gram_schmidt(char *matrix_filename) {
	printf("Test orthonormalization\n");
	int size[2];
	double *mat;
	read_matrix(matrix_filename, size, &mat);
	double (*A)[size[1]] = (double (*) []) mat;

	orthonormalize(size[0], size[1], A);

	for (int i = 0; i<size[0]; i++) {

		for ( int j = 0; j<size[1]; j++) {
			printf("[%f]", A[i][j]);
		}
		printf("\n");
	}

	printf("\n");
	for (int i = 0; i<3; i++) {
		double tmp2 = 0;
		double tmp = 0;
		for ( int j = 0; j<3; j++) {
			tmp += pow(A[j][i], 2);
			tmp2 += A[j][i]*A[j][(i + 1)%3];
		}
		if ((float) tmp != 1.0)
			printf("Wrong normalization! Expected 1 got : %lf\n", tmp);
		if (abs(tmp2) > 1e-14)
			printf("Wrong orthogonality! Expected 0 got : %lf\n", tmp2);
	}
	free(mat);
}

void test_read() {
	int size[2];
	int M = 5;
	double *mat;
	read_mtx("../test/matrices/bcspwr01.mtx", size, &mat);
	double (*A)[size[1]] = (double (*) []) mat;
	int N = size[0];

	double (*q)[M] = (double (*)[M]) malloc(sizeof(double)*N*M);
	init_q(N, M, q);

	mis(N, M, A, q, 150);

	double *expected_evec;
	int size_evec[2];
	read_matrix("../test/matrices/bcspwr01.eve", size_evec, &expected_evec);
	double (*expected)[size_evec[1]] = (double (*) []) expected_evec;

	for (int i = 0; i<size_evec[0]; i++) {
		double *vector = A[i];
		double *result = expected[i];
		double diff = scalar_product(N, vector, result)/(scalar_product(N, result, result)*scalar_product(N, vector, vector)); 
		printf("Distance : %f\n", diff);

	}
}

int main() {
	test_read();
	return 0;
}

