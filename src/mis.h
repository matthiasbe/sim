void orthonormalize(int N, int M, double A[M][N]);
void init(int N, int M, double (*A)[N], double (*q)[N]);
void mis(int N, int M, double (*A)[N], double (*q)[M], int iter);
double scalar_product(int N, double u[N], double v[N]);
void matrix_product(int M, int N, int P, double A[M][N], double B[N][P], double C[M][P]);
void transpose(int M, int N, double A[M][N], double T[N][M]);
void projection(int N, int M, double A[N][N], double Z[N][M], double B[M][M]);

