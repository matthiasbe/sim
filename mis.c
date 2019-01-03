/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>

double scalar_product(double u[3], double v[3]) {
	double res = 0;

	for ( int i; i<3; i++) {
		res += u[i]*v[i];
		printf("res = %f\n", res);
	}

printf("resres = %f\n", res);
	return res;
}

void orthonormalize(double A[3][3]) {
	
	double temp[3];
	double norm;

	for(int i = 0; i<3; i++) {
		printf("Traitement du vecteur numéro %d\n", i);
		for (int j = 0; j<3; j++) {
			temp[j] = A[i][j];
		}

		for (int k = 0; k<i; k++) {

			printf("calcul de la project sur l'axe %d\n", k);

			norm = scalar_product(A[i], A[k]);
			printf("%f\n", norm);
			norm /= scalar_product(A[k], A[k]);

			printf("%f\n", norm);

			for (int j = 0; j<3; j++) {
				temp[j] -= norm * A[i][j];
			}

		}

		for (int j = 0; j<3; j++) {
			A[i][j] = temp[j];
		}
	}
}

int main() {

        double A[3][3] = {
                {1,2,3},
                {2,5,3},
                {3,3,15}
        };

	// Initialise S0
        double q[2][3] = {
		{4.3,7.87,4.64},
		{4.3,-7.87,-4.64}
	};

	// Temp vector
        double v[2][3];
        double alpha[2];

	// A^k*v serie calculation
        for(int n = 0; n < 20; n++) {

		alpha[0] = 0;
		alpha[1] = 0;

		// Compute normalization dividende
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i<3;i++) {
				if (alpha[i] < abs(q[i][j]))
					alpha[i] = q[i][j];
			}
                }

		// v = A * Q
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
		orthonormalize(A);
		
        }

}


