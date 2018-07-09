#pragma once
#include <string>
double* Smoothing_Spine(double* X, double* F_wrong, double* P, const int& N);
double Cubic_Spline(double x, double* X, double* F, double* M, const int& N);
double* Cubic_Spline_Vector(double* X, double* F, const int& N);
double* approximation_coefs(double* xSet, double* fSet, const int& setSize, const int& polinomialPower = 5, std::string functions = "polinome", double SystemTolerance = 0.01);
double* approximation_Chebishov(double* xSet, double* fSet, const int& setSize, const int& polinomialPower = 4);
int approximation_optimal_N(double* xSet, double* FSet, const int& setSize, std::string approx = "polinome", const double& eps = 0.01);
double interpolation_Lagrange(double x, double* xSet, double* fSet, const int& setSize);
bool PolinomsComparator(double c0, double c1, double c2, double c3, double c4,
	double d0, double d1, double d2, double d3, double d4,
	double xFrom, double xTo);
double** create_approx_polinomes(double **matrix, const int& hight, const int& width, const int& polinome_power, std::string approximationFunc = "polinome");