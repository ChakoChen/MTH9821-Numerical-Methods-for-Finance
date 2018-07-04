//BarrierOption
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "BarrierOption.hpp"
using namespace std;

//To do Tridiagonal LU decomposition on A;
LUResult Tridiagonal_LU(vector<vector<double>> A) {
	int sz = A.size();
	vector<vector<double>> L(sz, vector<double>(sz, 0));
	vector<vector<double>> U(sz, vector<double>(sz, 0));

	for (int i = 0; i < sz - 1; i++) {
		L[i][i] = 1;
		L[i + 1][i] = A[i + 1][i] / A[i][i];
		U[i][i] = A[i][i];
		U[i][i + 1] = A[i][i + 1];
		A[i + 1][i + 1] = A[i + 1][i + 1] - L[i + 1][i] * U[i][i + 1];
	}
	L[sz - 1][sz - 1] = 1;
	U[sz - 1][sz - 1] = A[sz - 1][sz - 1];

	LUResult Result = { L, U };
	return Result;
}

//pdf of standard normal
double pdf(double x) {
	return exp(-pow(x, 2) / 2) / (sqrt(2 * M_PI));
}

//cdf of standard normal
double cdf(double x) {
	return  0.5 + erf(x / sqrt(2)) / 2;
}

//Domain Discretization
vector<double> BarrierOption::DomainDiscret(double alpha, int M) {
	double x_compute = log(S / K);
	double x_left = log(B / K);
	double tao_final = 0.5*T*pow(sigma, 2);
	double delta_tao = tao_final / M;

	double delta_x_temp = sqrt(delta_tao / alpha);
	int N_left = floor((x_compute - x_left) / delta_x_temp);
	double delta_x = (x_compute - x_left) / N_left;

	alpha = delta_tao / (pow(delta_x, 2));

	double x_right = log(S / K) + (r - q - 0.5*pow(sigma, 2))*T + 3 * sigma*sqrt(T);
	int N_right = ceil((x_right - x_compute) / delta_x);
	int N = N_left + N_right;
	x_right = x_compute + N_right*delta_x;

	//alpha, x_left, x_right, N, delta_x, delta_tao
	vector<double> Result = { alpha, x_left, x_right, double(N), delta_x, delta_tao};
	return Result;
}

//Black_Scholes for Call
vector<double> BarrierOption::Black_Scholes_Vanilla(string type) {
	double d1 = (log(S / K) + (r - q + pow(sigma, 2)*0.5)*T) / (sigma*sqrt(T));
	double d2 = d1 - sigma*sqrt(T);

	double calltheta = q*S*exp(-q*T)*cdf(d1) - r*K*exp(-r*T)*cdf(d2) -
		K*exp(-r*T)*pdf(d2)*sigma / (2 * sqrt(T));

	vector<double> Result;
	if (type.compare("PUT") == 0) {
		double P = K*exp(-r*T)*cdf(-d2) - S*exp(-q*T)*cdf(-d1);			//PValue
		double DeltaP = -exp(-q*T)*cdf(-d1);							//PDelta
		double GammaP = exp(-q*T)*pdf(d1) / (S*sigma*sqrt(T));			//PGamma
		double ThetaP = calltheta - q*S*exp(-q*T) + r*K*exp(-r*T);		//PTheta
		double VegaP = S*exp(-q*T)*pdf(d1)*sqrt(T);
		Result = { P, DeltaP, GammaP, ThetaP, VegaP };
	}
	else if (type.compare("CALL") == 0) {
		double C = S*exp(-q*T)*cdf(d1) - K*exp(-r*T)*cdf(d2);			//Cvalue
		double DeltaC = exp(-q*T)*cdf(d1);								//CDelta
		double GammaC = exp(-q*T)*pdf(d1) / (S*sigma*sqrt(T));			//CGamma
		double ThetaC = calltheta;										//CTheta
		double VegaC = S*exp(-q*T)*pdf(d1)*sqrt(T);
		Result = { C, DeltaC, GammaC, ThetaC, VegaC };
	}
	return Result;
}

//Black_Scholes Value
vector<double> BarrierOption::Black_Scholes_Barrier() {
	double a = (r - q) / pow(sigma, 2) - 0.5;
	double C1 = Black_Scholes_Vanilla("CALL")[0];
	BarrierOption Opt(pow(B, 2) / S, K, B, T, q, r, sigma);
	double C2 = Opt.Black_Scholes_Vanilla("CALL")[0];
	vector<double> Result = { C1 - pow(B / S, 2 * a)*C2 };
	return Result;

}

//Forward Euler
vector<double> BarrierOption::ForwardEuler(double alpha, int M, string Umatrix) {
	//Domain Discretization
	double x_compute = log(S / K);
	double x_left = log(B / K);
	double tao_final = 0.5*T*pow(sigma, 2);
	double delta_tao = tao_final / M;
	double delta_x_temp = sqrt(delta_tao / alpha);
	int N_left = floor((x_compute - x_left) / delta_x_temp);
	double delta_x = (x_compute - x_left) / N_left;
	alpha = delta_tao / (pow(delta_x, 2));
	double x_right = log(S / K) + (r - q - 0.5*pow(sigma, 2))*T + 3 * sigma*sqrt(T);
	int N_right = ceil((x_right - x_compute) / delta_x);
	int N = N_left + N_right;
	x_right = x_compute + N_right*delta_x;

	//change of variables:
	double a = (r - q) / pow(sigma, 2) - 1.0 / 2;
	double b = pow((r - q) / pow(sigma, 2) + 0.5, 2) + 2 * q / pow(sigma, 2);

	//construct the result for each point
	vector<vector<double>> Result(M + 1, vector<double>(N + 1));

	//calculte for tao = 0 (m=0)
	for (int i = 0; i < N + 1; i++) {
		double x_i = x_left + i*delta_x;
		Result[0][i] = K*exp(a*x_i)*max(exp(x_i) - 1, 0.0);
	}
	//calculate for x_left(n=0) and x_right(n=N)
	for (int i = 0; i < M + 1; i++) {
		double tao_i = i*delta_tao;
		//x_left
		Result[i][0] = 0;
		//x_right
		Result[i][N] = K*exp(a*x_right + b*tao_i)*(exp(x_right - 2 * q*tao_i / pow(sigma, 2)) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}

	//calculte for m = 1:M
	for (int m = 1; m <= M; m++) {
		for (int j = 1; j <= N - 1; j++) {
			Result[m][j] = alpha*Result[m - 1][j + 1] + (1 - 2 * alpha)*Result[m - 1][j] + alpha*Result[m - 1][j - 1];
		}
	}
	
	//congregate result
	double u_value = Result[M][N_left];
	double Option_Value = u_value*exp(-a*x_compute - b*tao_final);
	double error_point = abs(Option_Value - Black_Scholes_Barrier()[0]);

	double x_1 = x_compute - delta_x;
	double x0 = x_compute;
	double x1 = x_compute + delta_x;

	double S_1 = K*exp(x_1);
	double S0 = K*exp(x0);
	double S1 = K*exp(x1);

	double V_1 = exp(-a*x_1 - b*tao_final)*Result[M][N_left - 1];
	double V0 = exp(-a*x0 - b*tao_final)*Result[M][N_left];
	double V1 = exp(-a*x1 - b*tao_final)*Result[M][N_left + 1];
	
	double Delta = (V1 - V_1) / (S1 - S_1);
	double Gamma = 2 * ((S0 - S_1)*V1 - (S1 - S_1)*V0 + (S1 - S0)*V_1) / ((S0 - S_1)*(S1 - S0)*(S1 - S_1));
	double delta_t = 2 * delta_tao / pow(sigma, 2);
	double V_appro = exp(-a*x0 - b*(tao_final - delta_tao))*Result[M - 1][N_left];
	double Theta = (Option_Value - V_appro) / delta_t;

	if (Umatrix.compare("Yes") == 0) {
		for (int m = 0; m <= M; m++) {
			for_each(Result[m].begin(), Result[m].end(), [](double i) { cout << i << ", " ; });
			cout << endl;
		}
	}
	return { u_value, Option_Value, error_point, Delta, Gamma, Theta };
}

//Backward Euler
vector<double> BarrierOption::BackwardEuler_LU(double alpha, int M, string Umatrix) {
	//Domain Discretization
	double x_compute = log(S / K);
	double x_left = log(B / K);
	double tao_final = 0.5*T*pow(sigma, 2);
	double delta_tao = tao_final / M;
	double delta_x_temp = sqrt(delta_tao / alpha);
	int N_left = floor((x_compute - x_left) / delta_x_temp);
	double delta_x = (x_compute - x_left) / N_left;
	alpha = delta_tao / (pow(delta_x, 2));
	double x_right = log(S / K) + (r - q - 0.5*pow(sigma, 2))*T + 3 * sigma*sqrt(T);
	int N_right = ceil((x_right - x_compute) / delta_x);
	int N = N_left + N_right;
	x_right = x_compute + N_right*delta_x;

	//change of variables:
	double a = (r - q) / pow(sigma, 2) - 0.5;
	double b = pow((r - q) / pow(sigma, 2) + 0.5, 2) + 2 * q / pow(sigma, 2);

	//construct maxtrix A(size: N-1)
	vector<vector<double>>  A(N - 1, vector<double>(N - 1, 0));
	for (int i = 1; i < N - 2; i++) {
		A[i][i] = 1 + 2 * alpha;
		A[i][i - 1] = -1 * alpha;
		A[i][i + 1] = -1 * alpha;
	}
	A[0][0] = 1 + 2 * alpha;
	A[0][1] = -1 * alpha;
	A[N - 2][N - 2] = 1 + 2 * alpha;
	A[N - 2][N - 3] = -1 * alpha;

	//Get Tridiagonal_LU decomposition of A
	LUResult LandU = Tridiagonal_LU(A);
	vector<vector<double>> L = LandU.Lmatrix;
	vector<vector<double>> U = LandU.Umatrix;

	//construct the result for each point
	vector<vector<double>> Result(M + 1, vector<double>(N + 1));

	//calculte for tao = 0 (m=0)
	for (int i = 0; i < N + 1; i++) {
		double x_i = x_left + i*delta_x;
		Result[0][i] = K*exp(a*x_i)*max(exp(x_i) - 1, 0.0);
	}
	//calculate for x_left(n=0) and x_right(n=N)
	for (int i = 0; i < M + 1; i++) {
		double tao_i = i*delta_tao;
		//x_left
		Result[i][0] = 0;
		//x_right
		Result[i][N] = K*exp(a*x_right + b*tao_i)*(exp(x_right - 2 * q*tao_i / pow(sigma, 2)) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}

	//calculate for m=1:M
	for (int m = 1; m < M + 1; m++) {

		//construct vector b(size:N-1);
		vector<double> b(N - 1);
		for (int i = 1; i< N - 2; i++) {
			b[i] = Result[m - 1][i + 1];
		}
		b[0] = Result[m - 1][1] + alpha*Result[m][0];
		b[N - 2] = Result[m - 1][N - 1] + alpha*Result[m][N];

		//Solve for LUx=b : find y s.t. Ly=b (size of y: N-1)
		vector<double> Y(N - 1);
		Y[0] = b[0];
		for (int i = 1; i < N - 1; i++) {
			Y[i] = b[i] - L[i][i - 1] * Y[i - 1];
		}

		//Solve for Ux=y £º find x s.t. Ux= y(size of x: N-1)
		Result[m][N - 1] = Y[N - 2] / U[N - 2][N - 2];
		for (int j = N - 3; j >= 0; j--) {
			Result[m][j + 1] = (Y[j] - U[j][j + 1] * Result[m][j + 2]) / U[j][j];
		}
	}

	//congregate result
	double u_value = Result[M][N_left];
	double Option_Value = u_value*exp(-a*x_compute - b*tao_final);
	double error_point = abs(Option_Value - Black_Scholes_Barrier()[0]);

	double x_1 = x_compute - delta_x;
	double x0 = x_compute;
	double x1 = x_compute + delta_x;

	double S_1 = K*exp(x_1);
	double S0 = K*exp(x0);
	double S1 = K*exp(x1);

	double V_1 = exp(-a*x_1 - b*tao_final)*Result[M][N_left - 1];
	double V0 = exp(-a*x0 - b*tao_final)*Result[M][N_left];
	double V1 = exp(-a*x1 - b*tao_final)*Result[M][N_left + 1];

	double Delta = (V1 - V_1) / (S1 - S_1);
	double Gamma = 2 * ((S0 - S_1)*V1 - (S1 - S_1)*V0 + (S1 - S0)*V_1) / ((S0 - S_1)*(S1 - S0)*(S1 - S_1));
	double delta_t = 2 * delta_tao / pow(sigma, 2);
	double V_appro = exp(-a*x0 - b*(tao_final - delta_tao))*Result[M - 1][N_left];
	double Theta = (Option_Value - V_appro) / delta_t;

	if (Umatrix.compare("Yes") == 0) {
		for (int m = 0; m <= M; m++) {
			for_each(Result[m].begin(), Result[m].end(), [](double i) { cout << i << ", " ; });
			cout << endl;
		}
	}
	
	return { u_value, Option_Value, error_point, Delta, Gamma, Theta };

}

//Crank Nicolson
vector<double> BarrierOption::CrankNicolson_SOR(double alpha, int M, double w, double tol) {
	//Domain Discretization
	double x_compute = log(S / K);
	double x_left = log(B / K);
	double tao_final = 0.5*T*pow(sigma, 2);
	double delta_tao = tao_final / M;
	double delta_x_temp = sqrt(delta_tao / alpha);
	int N_left = floor((x_compute - x_left) / delta_x_temp);
	double delta_x = (x_compute - x_left) / N_left;
	alpha = delta_tao / (pow(delta_x, 2));
	double x_right = log(S / K) + (r - q - 0.5*pow(sigma, 2))*T + 3 * sigma*sqrt(T);
	int N_right = ceil((x_right - x_compute) / delta_x);
	int N = N_left + N_right;
	x_right = x_compute + N_right*delta_x;

	//change of variables:
	double a = (r - q) / pow(sigma, 2) - 0.5;
	double b = pow((r - q) / pow(sigma, 2) + 0.5, 2) + 2 * q / pow(sigma, 2);

	//construct maxtrix A(size: N-1)
	vector<vector<double>>  A(N - 1, vector<double>(N - 1, 0));
	for (int i = 1; i < N - 2; i++) {
		A[i][i] = 1 + alpha;
		A[i][i - 1] = -0.5 * alpha;
		A[i][i + 1] = -0.5 * alpha;
	}
	A[0][0] = 1 + alpha;
	A[0][1] = -0.5 * alpha;
	A[N - 2][N - 2] = 1 + alpha;
	A[N - 2][N - 3] = -0.5 * alpha;

	//construct the result for each point
	vector<vector<double>> Result(M + 1, vector<double>(N + 1));

	//calculte for tao = 0 (m=0)
	for (int i = 0; i < N + 1; i++) {
		double x_i = x_left + i*delta_x;
		Result[0][i] = K*exp(a*x_i)*max(exp(x_i) - 1, 0.0);
	}
	//calculate for x_left(n=0) and x_right(n=N)
	for (int i = 0; i < M + 1; i++) {
		double tao_i = i*delta_tao;
		//x_left
		Result[i][0] = 0;
		//x_right
		Result[i][N] = K*exp(a*x_right + b*tao_i)*(exp(x_right - 2 * q*tao_i / pow(sigma, 2)) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}

	//calculate for m=1:M
	for (int m = 1; m <= M; m++) {
		//construct vector b(size:N-1);
		vector<double> b(N - 1);
		for (int i = 1; i <= N - 3; i++) {
			b[i] = 0.5*alpha*Result[m - 1][i + 2] + (1 - alpha)*Result[m - 1][i + 1] + 0.5*alpha*Result[m - 1][i];
		}
		b[0] = 0.5*alpha*Result[m - 1][2] + (1 - alpha)*Result[m - 1][1]
			+ 0.5*alpha*Result[m - 1][0] + 0.5*alpha*Result[m][0];
		b[N - 2] = 0.5*alpha*Result[m - 1][N] + (1 - alpha)*Result[m - 1][N - 1]
			+ 0.5*alpha*Result[m - 1][N - 2] + 0.5*alpha*Result[m][N];

		//Initianization for X_New
		vector<double> X_Old(N + 1);
		vector<double> X_New(N + 1);
		//PAY ATTENTTION TO BOUNDARY CONDITIONS WHEN USING SOR
		X_Old[0] = 0;
		X_Old[N] = 0;
		X_New[0] = 0;
		X_New[N] = 0;
		for (int i = 1; i <= N - 1; i++) {
			X_Old[i] = S - K;
		}

		//Using SOR to find Result[m][]
		double diff = tol + 1;		//to ensure that it will enter the iteration
		while (diff > tol) {
			//Calculate New X
			for (int j = 1; j <= N - 1; j++) {
				X_New[j] = (1 - w)*X_Old[j] + 0.5*w*alpha*(X_New[j - 1] + X_Old[j + 1]) / (1 + alpha) + w*b[j - 1] / (1 + alpha);
			}
			//calculate diff 
			diff = 0;
			for (int i = 1; i <= N - 1; i++) {
				diff += pow(X_New[i] - X_Old[i], 2);
			}
			diff = sqrt(diff);
			//Switch
			X_Old = X_New;
		}

		//assign value of X_New to Result[m][]
		for (int j = 1; j <= N-1; j++) {
			Result[m][j] = X_Old[j];
		}
	}

	//congregate result
	double u_value = Result[M][N_left];
	double Option_Value = u_value*exp(-a*x_compute - b*tao_final);
	double error_point = abs(Option_Value - Black_Scholes_Barrier()[0]);

	double x_1 = x_compute - delta_x;
	double x0 = x_compute;
	double x1 = x_compute + delta_x;

	double S_1 = K*exp(x_1);
	double S0 = K*exp(x0);
	double S1 = K*exp(x1);

	double V_1 = exp(-a*x_1 - b*tao_final)*Result[M][N_left - 1];
	double V0 = exp(-a*x0 - b*tao_final)*Result[M][N_left];
	double V1 = exp(-a*x1 - b*tao_final)*Result[M][N_left + 1];

	double Delta = (V1 - V_1) / (S1 - S_1);
	double Gamma = 2 * ((S0 - S_1)*V1 - (S1 - S_1)*V0 + (S1 - S0)*V_1) / ((S0 - S_1)*(S1 - S0)*(S1 - S_1));
	double delta_t = 2 * delta_tao / pow(sigma, 2);
	double V_appro = exp(-a*x0 - b*(tao_final - delta_tao))*Result[M - 1][N_left];
	double Theta = (Option_Value - V_appro) / delta_t;

	return { u_value, Option_Value, error_point, Delta, Gamma, Theta };



}
