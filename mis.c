/************************************
 * Méthode des itérations simultanées
 ************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3
#define M 2

double scalar_product(double u[N], double v[N]) {
	double res = 0;

	for ( int i = 0; i<N; i++) {
		res += u[i]*v[i];
	}

	return res;
}

void orthonormalize(double A[M][N]) {
	
	double temp[N];
	double norm, temp_norm;

	for(int i = 0; i<M; i++) {
//		printf("Traitement du vecteur numéro %d\n", i);
		for (int j = 0; j<N; j++) {
			temp[j] = A[i][j];
		}

		for (int k = 0; k<i; k++) {

//			printf("calcul de la project sur l'axe %d\n", k);

			norm = scalar_product(A[i], A[k]);
			printf("(u,v) = %f\n", norm);
			norm /= scalar_product(A[k], A[k]);

			printf("(u,v) / (u,u) = %f\n", norm);

			for (int j = 0; j<N; j++) {
				temp[j] -= norm * A[i][j];
			}
		}
		for (int j = 0; j<N; j++) {
			A[i][j] = temp[j];
		}
	}

	for(int i = 0; i<M; i++) {
		temp_norm = 0;
		for ( int k = 0; k<N; k++) {
			temp_norm += temp[k] * temp[k];
		}

		temp_norm = sqrt(temp_norm);

		for (int j = 0; j<N; j++) {
			A[i][j] /= temp_norm;
		}
	}


}

int main() {

        double A[N][N] = {
                {1,2,3},
                {2,5,3},
                {3,3,15}
        };

	// Initialise S0
        double q[M][N] = {
		{0,3,1},
		{5,-2,3}
	};

	// Temp vector
        double v[M][N];

	// A^k*v serie calculation
        for(int n = 0; n < 5; n++) {

		// v = A * Q
		for (int k = 0; k<M; k++) {
			for (int i = 0; i<N;i++) {
				v[k][i] = 0.0;
				for (int j = 0; j<3;j++) {
					v[k][i] += q[k][j]*A[i][j];
				}
				v[k][i] = v[k][i];
			}
		}
		// q = v
                for(int i = 0; i<M; i++) {
			for(int j = 0; j<N; j++) {
				printf("[%f]",v[i][j]);
			}
			printf("\n");
                }
                printf("\n");

		printf("Gram-Schmidt ...\n");
		orthonormalize(v);

		// q = v
                for(int i = 0; i<M; i++) {
			for(int j = 0; j<N; j++) {
				q[i][j] = v[i][j];
				printf("[%f]",q[i][j]);
			}
			printf("\n");
                }
                printf("\n");
		
        }

}


