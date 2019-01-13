#include <stdlib.h>
#include <stdio.h>
#include "mis.h"

void init_q(int N, int M, double q[M][N]){
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j<N; j++) {
			q[i][j] = (double) rand();
		}
	}
	
}

void test3x3() {
	double A[3][3] = {
		{154.61455, 123.1110, 231564.1},
		{231564.1 , 11.0005 , 789.12  },
		{231564.1 , 789.12  , 0.00254 }
	};

//	double expected_eval[3] = {232097.75470006, -231154.01191683, -778.12519323};
//
//	double expected_evec[3][3] = {
//			{0.576916252, 0.577988900, 0.0034060215},
//			{0.577580793, -0.577016915, 0.0999994065},
//			{0.577553517, -0.577044462, 0.000517927165}
//	};

	double computed_evec[3][2];
	init_q(3, 2, computed_evec);

	mis(3, 2, A, computed_evec, 5);


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
