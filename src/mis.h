#include <mpi.h>

void orthonormalize(int N, int M, double A[M][N]);
void init(int N, int M, double (*A)[N], double (*q)[N]);
void mis(int N, int M, double (*A)[N], double (*q)[M], int iter, int precision, MPI_Comm comm);
double scalar_product(int N, double u[N], double v[N]);
void matrix_product(int M, int N, int P, double A[M][N], double B[N][P], double C[M][P], MPI_Comm comm);
void transpose(int M, int N, double A[M][N], double T[N][M]);
void projection(int N, int M, double A[N][N], double Z[N][M], double B[M][M], MPI_Comm comm);
void print_matrix(int N, int M, double A[N][M]);
void init_q(int N, int M, double q[M][N]);
void compute_submatrix(int psize[2], int rank, int N, int M, int min[2], int max[2], MPI_Comm comm);
