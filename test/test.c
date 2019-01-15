#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../src/mis.h"
#include "../src/read_file.h"

void init_q(int N, int M, double q[M][N]){
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j<N; j++) {
			q[i][j] = (double) rand();
		}
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

	mis(3, 3, (double (*)[]) A, computed_evec, 50);

	for (int i = 0; i<3; i++) {

		for ( int j = 0; j<3; j++) {
			printf("[%f]", computed_evec[i][j]);
		}
		printf("\n");
	}
}

void test_gram_schmidt(char *matrix_filename) {
	printf("Test orthonormalization\n");
	int size[2];
	double *mat;
	read_matrix(matrix_filename, size, &mat);
	double (*A)[size[1]] = (double (*) []) mat;

	orthonormalize(size[0], size[1], A);

	for (int i = 0; i<3; i++) {

		for ( int j = 0; j<3; j++) {
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
			tmp += A[j][i]*A[j][(i + 1)%3];
		}
		if ((float) tmp != 1.0)
			printf("Expected 1 got : %f\n", tmp, 30);
		if (tmp2 != 0.0)
			printf("Expected 0 got : %f\n", tmp2, 30);
	}
}

int main() {
	//test3x3();
	test_gram_schmidt("../test/matrices/simple3x3");
	test_gram_schmidt("../test/matrices/3x3");
	return 0;
}

