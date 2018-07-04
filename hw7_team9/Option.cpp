#define _USE_MATH_DEFINES
#include <math.h>
#include<vector>
#include<iostream>
#include<algorithm>
#include"Option.hpp"
using namespace std;

//pdf of standard normal
double pdf(double x) {
	return exp(-pow(x, 2) / 2) / (sqrt(2 * M_PI));
}

//cdf of standard normal
double cdf(double x) {
	return  0.5 + erf(x / sqrt(2)) / 2;
}

//To do Tridiagonal LU decomposition on A;
LUResult Tridiagonal_LU(vector<vector<double>> A) {
	int sz = A.size();
	vector<vector<double>> L(sz, vector<double>(sz, 0));
	vector<vector<double>> U(sz, vector<double>(sz, 0));

	for (int i = 0; i < sz-1; i++) {
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

//return be values calculated by Black_Scholes
vector<double> Option::BlackScholes() {
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

//return be values calculated by Backward Euler Method
PriceResult Option::BackWardEuler(int M, double alpha_temp) {
	
	//change of variables
	double tao_final = T*pow(sigma, 2) / 2;
	double x_left = log(S / K) + (r - q - pow(sigma, 2) / 2)*T - 3 * sigma*sqrt(T);
	double x_right = log(S / K) + (r - q - pow(sigma, 2) / 2)*T + 3 * sigma*sqrt(T);
	double delta_Tao = tao_final / M;
	int N = floor((x_right - x_left) / sqrt(delta_Tao / alpha_temp));
	double delta_x = (x_right - x_left) / N;
	double alpha = delta_Tao / pow(delta_x, 2);

	double a = (r - q) / pow(sigma, 2) - 1.0 / 2;
	double b = pow((r - q) / pow(sigma, 2) + 1.0 / 2, 2) + 2 * q / pow(sigma, 2);


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
		Result[0][i] = K*exp(a*x_i)*max(1 - exp(x_i), 0.0);
	}
	//calculate for x_left(n=0) and x_right(n=N)
	for (int i = 0; i < M + 1; i++) {
		double tao_i = i*delta_Tao;
		//x_left
		Result[i][0] = K*exp(a*x_left + b*tao_i)*(exp(-2 * tao_i*r / pow(sigma, 2)) - exp(x_left - 2 * q*tao_i / pow(sigma, 2)));
		//x_tight
		Result[i][N] = 0;
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
		Result[m][N - 1] = Y[N - 2] / U[N -2][N - 2];
		for (int j = N - 3; j >= 0; j--) {
			Result[m][j+1] = (Y[j] - U[j][j + 1] * Result[m][j + 2]) / U[j][j];
		}
	}
	
	//Analysis Convergence
	double x_comp = log(S / K);
	int i = floor((x_comp - x_left) / delta_x);
	double x_1 = i*delta_x + x_left;
	double x_2 = (i + 1)*delta_x + x_left;

	double S_1 = K*exp(x_1);
	double S_2 = K*exp(x_2);
	double V_1 = exp(-a*x_1 - b*tao_final)*Result[M][i];
	double V_2 = exp(-a*x_2 - b*tao_final)*Result[M][i + 1];
	double V_appr1 = ((S_2 - S)*V_1 + (S - S_1)*V_2) / (S_2 - S_1);
	double V_bs = BlackScholes()[0];
	double error_pointwise = abs(V_appr1 - V_bs);

	double u_com = ((x_2 - x_comp)*Result[M][i] + (x_comp - x_1)*Result[M][i + 1]) / (x_2 - x_1);
	double V_appr2 = exp(-a*x_comp - b*tao_final)*u_com;
	double error_pointwise2 = abs(V_appr2 - V_bs);

	//RMS
	int count_num = 0;
	double count_sum = 0;
	for (int k = 0; k <= N; k++) {
		double xk = x_left + k*delta_x;
		double Sk = K*exp(xk);
		Option tempOp(Sk, K, T, r, q, sigma, type);
		double tempV_bs = tempOp.BlackScholes()[0];

		if (tempV_bs/S > 0.00001) {
			count_num += 1;
			double tempV_appr = exp(-a*xk - b*tao_final)*Result[M][k];
			count_sum += pow(tempV_appr - tempV_bs, 2) / pow(tempV_bs, 2);
		}
	}
	double error_RMS = sqrt(count_sum / count_num);

	//Greeks -- Delta
	double delta = (V_2 - V_1) / (S_2 - S_1);
	//Greeks -- Gamma
	double x_0 = (i - 1)*delta_x + x_left;
	double x_3 = (i + 2)*delta_x + x_left;
	double S_0 = K*exp(x_0);
	double S_3 = K*exp(x_3);
	double V_0 = exp(-a*x_0 - b*tao_final)*Result[M][i - 1];
	double V_3 = exp(-a*x_3 - b*tao_final)*Result[M][i];
	double gamma = ((V_3 - V_2) / (S_3 - S_2) - (V_1 - V_0) / (S_1 - S_0)) / ((S_3 + S_2) / 2 - (S_1 + S_0) / 2);
	//Greeks -- Theta
	double delta_t = 2 * delta_Tao / pow(sigma, 2);
	double V_1t = exp(-a*x_1 - b*(tao_final - delta_Tao))*Result[M - 1][i];
	double V_2t = exp(-a*x_2 - b*(tao_final - delta_Tao))*Result[M - 1][i + 1];
	double V_apprt = ((S_2 - S)*V_1t + (S - S_1)*V_2t) / (S_2 - S_1);
	double theta = (V_apprt - V_appr1) / delta_t;

	vector<double> Convergence = { error_pointwise, error_pointwise2, error_RMS, delta, gamma, theta };
	
	//Result is the 2-d vector contains values calculated for each nodes,
	//Convergence is the 1-d vector contains error_pointwise, error_pointwise2, error_RMS, Delta, Gamma, Theta
	PriceResult Final = { Result, Convergence };
	return Final;
}

//return be values calculated by CrankNicolson Method
PriceResult Option::CrankNicolson(int M, double alpha_temp) {

	//change of variables
	double tao_final = T*pow(sigma, 2) / 2;
	double x_left = log(S / K) + (r - q - pow(sigma, 2) / 2)*T - 3 * sigma*sqrt(T);
	double x_right = log(S / K) + (r - q - pow(sigma, 2) / 2)*T + 3 * sigma*sqrt(T);
	double delta_Tao = tao_final / M;
	int N = floor((x_right - x_left) / sqrt(delta_Tao / alpha_temp));
	double delta_x = (x_right - x_left) / N;
	double alpha = delta_Tao / pow(delta_x, 2);
	double a = (r - q) / pow(sigma, 2) - 1.0 / 2;
	double b = pow((r - q) / pow(sigma, 2) + 1.0 / 2, 2) + 2 * q / pow(sigma, 2);

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

	//Get Tridiagonal_LU decomposition of A
	LUResult LandU = Tridiagonal_LU(A);
	vector<vector<double>> L = LandU.Lmatrix;
	vector<vector<double>> U = LandU.Umatrix;

	//construct the result for each point
	vector<vector<double>> Result(M + 1, vector<double>(N + 1));

	//calculte for tao = 0 (m=0)
	for (int i = 0; i < N + 1; i++) {
		double x_i = x_left + i*delta_x;
		Result[0][i] = K*exp(a*x_i)*max(1 - exp(x_i), 0.0);
	}
	//calculate for x_left(n=0) and x_right(n=N)
	for (int i = 0; i < M + 1; i++) {
		double tao_i = i*delta_Tao;
		//x_left
		Result[i][0] = K*exp(a*x_left + b*tao_i)*(exp(-2 * tao_i*r / pow(sigma, 2)) - exp(x_left - 2 * q*tao_i / pow(sigma, 2)));
		//x_tight
		Result[i][N] = 0;
	}

	//calculate for m=1:M
	for (int m = 1; m < M + 1; m++) {
		//construct vector b(size:N-1);
		vector<double> b(N - 1);
		for (int i = 1; i < N - 2; i++) {
			b[i] = 0.5*alpha*Result[m - 1][i + 2] + (1 - alpha)*Result[m - 1][i + 1] + 0.5*alpha*Result[m - 1][i];
		}
		b[0] = 0.5*alpha*Result[m - 1][2] + (1 - alpha)*Result[m - 1][1] 
			+ 0.5*alpha*Result[m - 1][0] + 0.5*alpha*Result[m][0];
		b[N - 2] = 0.5*alpha*Result[m - 1][N] + (1 - alpha)*Result[m - 1][N - 1] 
			+ 0.5*alpha*Result[m - 1][N - 2] + 0.5*alpha*Result[m][N];

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
	//Analysis Convergence
	double x_comp = log(S / K);
	int i = floor((x_comp - x_left) / delta_x);
	double x_1 = i*delta_x + x_left;
	double x_2 = (i + 1)*delta_x + x_left;

	double S_1 = K*exp(x_1);
	double S_2 = K*exp(x_2);
	double V_1 = exp(-a*x_1 - b*tao_final)*Result[M][i];
	double V_2 = exp(-a*x_2 - b*tao_final)*Result[M][i + 1];
	double V_appr1 = ((S_2 - S)*V_1 + (S - S_1)*V_2) / (S_2 - S_1);
	double V_bs = BlackScholes()[0];
	double error_pointwise = abs(V_appr1 - V_bs);

	double u_com = ((x_2 - x_comp)*Result[M][i] + (x_comp - x_1)*Result[M][i + 1]) / (x_2 - x_1);
	double V_appr2 = exp(-a*x_comp - b*tao_final)*u_com;
	double error_pointwise2 = abs(V_appr2 - V_bs);

	//calculate RMS
	int count_num = 0;
	double count_sum = 0;
	for (int k = 0; k <= N; k++) {
		double xk = x_left + k*delta_x;
		double Sk = K*exp(xk);
		Option tempOp(Sk, K, T, r, q, sigma, type);
		double tempV_bs = tempOp.BlackScholes()[0];

		if (tempV_bs / S > 0.00001) {
			count_num += 1;
			double tempV_appr = exp(-a*xk - b*tao_final)*Result[M][k];
			count_sum += pow(tempV_appr - tempV_bs, 2) / pow(tempV_bs, 2);
		}
	}
	double error_RMS = sqrt(count_sum / count_num);

	//Greeks -- Delta
	double delta = (V_2 - V_1) / (S_2 - S_1);
	//Greeks -- Gamma
	double x_0 = (i - 1)*delta_x + x_left;
	double x_3 = (i + 2)*delta_x + x_left;
	double S_0 = K*exp(x_0);
	double S_3 = K*exp(x_3);
	double V_0 = exp(-a*x_0 - b*tao_final)*Result[M][i - 1];
	double V_3 = exp(-a*x_3 - b*tao_final)*Result[M][i];
	double gamma = ((V_3 - V_2) / (S_3 - S_2) - (V_1 - V_0) / (S_1 - S_0)) / ((S_3 + S_2) / 2 - (S_1 + S_0) / 2);
	//Greeks -- Theta
	double delta_t = 2 * delta_Tao / pow(sigma, 2);
	double V_1t = exp(-a*x_1 - b*(tao_final - delta_Tao))*Result[M - 1][i];
	double V_2t = exp(-a*x_2 - b*(tao_final - delta_Tao))*Result[M - 1][i + 1];
	double V_apprt = ((S_2 - S)*V_1t + (S - S_1)*V_2t) / (S_2 - S_1);
	double theta = (V_apprt - V_appr1) / delta_t;

	vector<double> Convergence = { error_pointwise, error_pointwise2, error_RMS, delta, gamma, theta };

	//Result is the 2-d vector contains values calculated for each nodes,
	//Convergence is the 1-d vector contains error_pointwise, error_pointwise2, error_RMS, Delta, Gamma, Theta
	PriceResult Final = { Result, Convergence };
	return Final;
}
