#include <stdexcept>
#include "LinearEquations.h"
#include "MatrixVectorMultiplication.h"
//typedef int var_type;	// basic sizes, can be changed
#define FLT_EPSILON 1.192093e-07
#define DBL_EPSILON 2.220446e-16

double* Cubic_Spline_Vector(double* X, double* F, const int& N) {
	// f : i = 0...N
	// X : i = 0...N
	double* M = new double[N];
	double* a = new double[N-2];
	double* b = new double[N-2];
	double* c = new double[N-2];
	double* d = new double[N-2];

	int i;
	for (i = 1; i < N - 1; ++i) {
		a[i-1] = (X[i] - X[i - 1]);
		b[i-1] = 2 * (X[i + 1] - X[i - 1]); //	2*(X[i] - X[i-1] + X[i+1] - X[i])
		c[i-1] = (X[i + 1] - X[i]);
		d[i-1] = 6 * ((F[i + 1] - F[i]) / (X[i + 1] - X[i]) - (F[i] - F[i - 1]) / (X[i] - X[i - 1]));
	}

	double* res = TridiagonalSolve(a, b, c, d, N - 2);
	M[0] = 0; M[N - 1] = 0;
	for (i = 1; i < N - 1; ++i) {
		M[i] = res[i];
	}
	delete[]a; delete[]b; delete[]c;
	delete[]d; delete[]res;
	return M;
}

double Cubic_Spline(double x, double* X, double* F, double* M, const int& N) {
	/*Returns an approximated value of f(x) = Q(x), x in [X[i], X[i-1]], 
												using natural cubic splines*/
	int i;
	for (i = 0; i < N - 1; ++i) {
		if (double(x - X[i]) >= 0 && double(X[i + 1] - x) > 0) break;
		else if (double(X[i + 1] - x) < DBL_EPSILON) {
			++i; break;
		}
	}

	double di, bi, ai, hi;
	hi = (X[i + 1] - X[i]);
	ai = F[i];
	di = (M[i] - M[i - 1]) / hi;
	bi = hi / 2 * M[i] - (hi*hi / 6)*di + (F[i] - F[i - 1]) / hi;
	double diff = (x - X[i]);
	return ai + bi * diff + M[i] / 2 * diff*diff + di / 6 * diff*diff*diff;
}

double* Smoothing_Spine(double* X, double* F_wrong, double* P, const int& N) {
	int i, j, k;
	double* h = new double[N];
	for (i = 1; i < N + 1; ++i) {
		h[i - 1] = X[i] - X[i - 1];
	}
	
	// H = (H1[i]	H2[i]	H3[i])
	// H has N - 1 rows and N + 1 collums
	double* H1 = new double[N];
	double* H2 = new double[N];
	double* H3 = new double[N];
	H1[0] = 1 / h[0]; H1[N - 1] = 1 / h[N - 2];
	H2[0] = (-1 / h[0] - 1 / h[1]);
	H2[N - 1] = (-1 / h[N - 2] - 1 / h[N - 1]);
	H3[0] = 1 / h[1]; H3[N - 1] = 1 / h[N - 1];

	// A = (A1[i]	A2[i]	A3[i])
	// A has N-1 rows and N-1 collums
	double* A1 = new double[N];
	double* A2 = new double[N];
	double* A3 = new double[N];
	A1[0] = 0;  A1[N - 1] = h[N - 1] / 3;
	A2[0] = (h[0] + h[1]) / 3;
	A2[N - 1] = (h[N - 2] + h[N - 1]) / 3;
	A3[0] = h[1] / 3;	 A3[N - 1] = 0;
	for (i = 1; i < N - 1; ++i) {
		H1[i] = 1 / h[i];
		H2[i] = (-1 / h[i] - 1 / h[i+1]);
		H3[i] = 1 / h[i + 1];

		A1[i] = h[i] / 3;
		A2[i] = (h[i] + h[i + 1]) / 3;
		A3[i] = h[i + 1] / 3;
	}

	delete[]h;

	double* h1 = new double[N];
	double* h2 = new double[N];
	double* h3 = new double[N];
	for (i = 0; i < N; ++i) {
		h1[i] = H1[i] / P[i];
		h2[i] = H2[i] / P[i + 1];
		h3[i] = H3[i] / P[i + 2];
	}

	double** resMatrix = new double*[N+2];
	double** HMatrix = new double*[N];
	for (i = 0; i < N + 2; ++i) {
		resMatrix[i] = new double[N];
	}
	for (i = 0; i < N; ++i) {
		HMatrix[i] = new double[N + 2];
	}

	// P^{-1} H*
	for (j = 0; j < N; ++j) {
		for (i = 0; i < j; ++i) {
			resMatrix[i][j] = 0;
		}
		resMatrix[j][j] = h1[j]; resMatrix[j+1][j] = h2[j]; resMatrix[j+2][j] = h3[j];
		for (i = j + 2; i < N + 2; ++i) {
			resMatrix[i][j] = 0;
		}
	}

	// H
	for (i = 0; i < N; ++i) {
		for (j = 0; j < i; ++j) {
			HMatrix[i][j] = 0;
		}
		HMatrix[i][i] = H1[i]; HMatrix[i][i+1] = H2[i]; HMatrix[i][i+2] = H3[i];
		for (j = i + 2; j < N + 2; ++j) {
			HMatrix[i][j] = 0;
		}
	}

	// (A + H P^{-1} H*)
	double** Result = Multiply_Matr(HMatrix, N, N + 2, resMatrix, N + 2, N);	//	N x N
	Result[0][0] += A2[0]; Result[0][1] += A3[0];
	for (i = 1; i < N - 1; ++i) {
		Result[i][i - 1] += A1[i];
		Result[i][i] += A2[i];
		Result[i + 1][i] += A3[i];
	}
	Result[N - 1][N - 2] = A1[N - 1]; Result[N - 1][N - 1] = A2[N - 1];

	//	H f_wr
	double* B = new double[N];
	for (i = 0; i < N; ++i) {
		B[i] = F_wrong[i] * H1[i] + F_wrong[i + 1] * H2[i] + F_wrong[i + 2] * H3[i];
	}

	delete[]h1;		delete[]h2;		delete[]h3;
	delete[]H1;		delete[]H2;		delete[]H3;
	delete[]A1;		delete[]A2;		delete[]A3;
	for (i = 0; i < N; ++i) {
		delete[]HMatrix[i];
	}
	delete[]HMatrix;

	double* m = GaussMethod(Result, B, N);
	double* nue = new double[N + 2];
	double* multiplied_vector = Multiply_vector(resMatrix, N + 2, N, m, N);
	double t = 0;
	for (i = 0; i < N + 2; ++i) {
		nue[i] = F_wrong[i] - multiplied_vector[i];
	}

	delete[]B;
	delete[]multiplied_vector;
	delete[]m;
	for (i = 0; i < N+2; ++i)
		delete[]resMatrix[i];

	delete[]resMatrix;
	for (i = 0; i < N; ++i) 
		delete[]Result[i];	
	delete[]Result;

	return nue;
}