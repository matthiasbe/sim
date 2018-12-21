/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>

int main() {

        double A[3][3] = {
                {1,2,3},
                {2,5,3},
                {3,3,0}
        };

	// Initialise S0
        double q[2][3] = {
		{1,0,0},
		{1,1,0}
	};

	// Temp vector
        double v[2][3];
        double alpha[2];

	// A^k*v serie calcultion
        for(int n = 0; n < 100; n++) {

		alpha[0] = 0;
		alpha[1] = 0;

		// Compute normalization dividende
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i<3;i++) {
				if (alpha[i] < q[i][j])
					alpha[i] = q[i][j];
			}
                }

		// v = A * q
		for (int k = 0; k<2; k++) {
			for (int i = 0; i<3;i++) {
				v[k][i] = 0.0;
				for (int j = 0; j<3;j++) {
					v[k][i] += q[k][j]*A[i][j];
				}
				v[k][i] = v[k][i] / alpha[k];
			}
		}

		// q = v
                for(int i = 0; i<2; i++) {
			for(int j = 0; j<3; j++) {
				q[i][j] = v[i][j];
				printf("[%f]",q[i][j]);
			}
			printf("\n");
                }
                printf("\n");
        }

}


