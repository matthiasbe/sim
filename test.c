#include <stdlib.h>
#include <stdio.h>
#include "mis.h"
#include "read_file.h"


void test3x3() {
	int size[2];
	double *mat;
	read_matrix("matrices/3x3", size, &mat);
	double (*A)[size[1]] = (double (*) []) mat;

	for ( int i = 0; i<3; i++) {

		for ( int j = 0; j<3; j++) {
			printf("[%f]", A[i][j]);
		}
		printf("\n");
	}
//	double expected_eval[3] = {232097.75470006, -231154.01191683, -778.12519323};
//
//	double expected_evec[3][3] = {
//			{0.576916252, 0.577988900, 0.0034060215},
//			{0.577580793, -0.577016915, 0.0999994065},
//			{0.577553517, -0.577044462, 0.000517927165}
//	};

	double computed_evec[3][2];

	mis(3, 2, (double (*)[]) A, computed_evec, 5);


	for ( int i = 0; i<3; i++) {

		for ( int j = 0; j<2; j++) {
			printf("[%f]", computed_evec[i][j]);
		}
		printf("\n");
	}
}

int main() {
	test3x3();
	return 0;
}
