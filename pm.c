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
        double q[3] = {1.0L,0.0L,0.0L};

	// Temp vector
        double v[3];
        double alpha;

	// A^k*v serie calcultion
        for(int k = 0; k < 100; k++) {

		// Compute normalization dividende
                for (int i = 0; i<3;i++) {
                        if (alpha < q[i])
                                alpha = q[i];
                }

		// Compute Sn = A * Sn-1
                for (int i = 0; i<3;i++) {
                        v[i] = 0.0;
                        for (int j = 0; j<3;j++) {
                                v[i] += q[j]*A[i][j];
                        }
                        v[i] = v[i] / alpha;
                }

		// Move v to q
                for(int i = 0; i<3; i++) {
                        q[i] = v[i];
                        printf("[%f]",q[i]);
                }
                printf("\n");
        }

}


