#include "FiniteDiff.hpp"
using namespace std;

tuple<vector<double>, vector<double>> FiniteDiff::DomainDiscrete_Div(double alpha_1, int M_1, double t_div, double alpha_2) {
	double S = getS();
	double K = getK();
	double T = getT();
	double sigma = getsigma();
	double q = getq();
	double r = getr();

	//Domain Discretization on [0, tao_div)
	double x_compute_ = log(S / K) + log(1 - q);
	double tao_final = 0.5*T*pow(sigma, 2);
	double tao_div = 0.5*(T - t_div)*pow(sigma, 2);
	double delta_tao1 = tao_div / M_1;
	double delta_x = sqrt(delta_tao1 / alpha_1);
	int N_left = ceil((x_compute_ - (log(S / K) + (r - 0.5*pow(sigma, 2))*T - 3 * sigma*sqrt(T)))/delta_x);
	int N_right = ceil(((log(S / K) + (r - 0.5*pow(sigma, 2))*T + 3 * sigma*sqrt(T)) - x_compute_) / delta_x);
	double x_left = x_compute_ - N_left*delta_x;
	double x_right = x_compute_ + N_right*delta_x;

	//Domain Discretization on [tao_div, tao_final]
	int M_2 = ceil((tao_final - tao_div) / (alpha_2*pow(delta_x, 2)));
	double delta_tao2 = (tao_final - tao_div) / M_2;
	alpha_2 = delta_tao2 / pow(delta_x, 2);


	vector<double> Result1 = { alpha_1, x_left, x_right, double(N_left), double(N_right), delta_x, delta_tao1 };
	vector<double> Result2 = { alpha_2, x_left - log(1 - q), x_right - log(1 - q), double(N_left), double(N_right), delta_x, delta_tao2, double(M_2)};
	return make_tuple(Result1, Result2);
}

//Forward Euler
tuple<vector<double>, mat, mat> FiniteDiff::ForwardEuler_Div(double alpha_1, int M_1, double t_div, double alpha_2) {
	double S = getS();
	double K = getK();
	double T = getT();
	double sigma = getsigma();
	double q = getq();
	double r = getr();

	//Domain Discretization
	auto res = DomainDiscrete_Div(alpha_1, M_1, t_div, alpha_2);
	alpha_1 = get<0>(res)[0];
	alpha_2 = get<1>(res)[0];
	double x_left = get<0>(res)[1];
	double x_right = get<0>(res)[2];
	int N_left = int(get<0>(res)[3]);
	int N = N_left + int(get<0>(res)[4]);
	double delta_x = get<0>(res)[5];
	double delta_tao_1 = get<0>(res)[6];
	double delta_tao_2 = get<1>(res)[6];
	int M_2 = int(get<1>(res)[7]);


	//change of variables:
	double a = r / pow(sigma, 2) - 0.5;
	double b = pow(a + 1, 2);
	double tao_div = 0.5*(T - t_div)*pow(sigma, 2);
	double tao_final = 0.5*T*pow(sigma, 2);

	//construct the result for each point
	mat Result_1(M_1 + 1, N + 1);
	Result_1.setZero();

	//calculte for tao = 0 (m=0)
	for (int i = 0; i < N + 1; i++) {
		double x_i = x_left + i*delta_x;
		Result_1(0, i) = K*exp(a*x_i)*max(exp(x_i) - 1, 0.0);
	}

	//calculate for x_left and x_right at [0, tao_div]
	for (int i = 0; i < M_1 + 1; i++) {
		double tao_i = i*delta_tao_1;
		//x_left
		Result_1(i, 0) = 0;
		//x_right
		Result_1(i, N) = K*exp(a*x_right + b*tao_i)*(exp(x_right) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}

	//calculte for m = 1:M_1 [0, tao_div]
	for (int m = 1; m <= M_1; m++) {
		for (int j = 1; j <= N - 1; j++) {
			Result_1(m, j) = alpha_1*Result_1(m - 1, j + 1) + (1 - 2 * alpha_1)*Result_1(m - 1, j) + alpha_1*Result_1(m - 1, j - 1);
		}
	}

	//construct the result for each point
	mat Result_2(M_2 + 1, N + 1);
	Result_2.setZero();

	//calculte for tao = tao_div (m=0)
	for (int i = 0; i < N + 1; i++) {
		Result_2(0, i) = pow(1 - q, -a)*Result_1(M_1, i);
	}

	//calculate for x_left and x_right at [tao_div, tao_final]
	for (int i = 1; i < M_2 + 1; i++) {
		double tao_i = tao_div + i*delta_tao_2;
		//x_left
		Result_2(i, 0) = 0;
		//x_right
		Result_2(i, N) = K*exp(a*(x_right - log(1 - q)) + b*tao_i)*(exp(x_right) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}


	//calculate for [tao_div, tao_final]
	for (int m = 1; m <= M_2; m++) {
		for (int j = 1; j <= N - 1; j++) {
			Result_2(m, j) = alpha_2*Result_2(m - 1, j + 1) + (1 - 2 * alpha_2)*Result_2(m - 1, j) + alpha_2*Result_2(m - 1, j - 1);
		}
	}

	//Pointwise convergence
	double x_compute = log(S / K);
	double x_1 = x_compute - delta_x;
	double x0 = x_compute;
	double x1 = x_compute + delta_x;

	double S_1 = K*exp(x_1);
	double S0 = K*exp(x0);
	double S1 = K*exp(x1);

	double V_1 = exp(-a*x_1 - b*tao_final)*Result_2(M_2, N_left - 1);
	double V0 = exp(-a*x0 - b*tao_final)*Result_2(M_2, N_left);
	double V1 = exp(-a*x1 - b*tao_final)*Result_2(M_2, N_left + 1);

	double Option_Value = V0;


	double Delta = (V1 - V_1) / (S1 - S_1);
	double Gamma = 2 * ((S0 - S_1)*V1 - (S1 - S_1)*V0 + (S1 - S0)*V_1) / ((S0 - S_1)*(S1 - S0)*(S1 - S_1));

	//Not sure to use delta_tao2 or not
	double delta_t = 2 * delta_tao_2 / pow(sigma, 2);
	double V_appro = exp(-a*x0 - b*(tao_final - delta_tao_2))*Result_2(M_2 - 1, N_left);
	double Theta = (Option_Value - V_appro) / delta_t;


	vector<double> VG = { Result_2(M_2, N_left), Option_Value, Delta, Gamma, Theta };
	return make_tuple(VG, Result_1, Result_2);
}

//CrankNicolson
tuple<vector<double>, mat, mat> FiniteDiff::CrankNicolson_Div(double alpha_1, int M_1, double t_div, double alpha_2) {
	double S = getS();
	double K = getK();
	double T = getT();
	double sigma = getsigma();
	double q = getq();
	double r = getr();

	//Domain Discretization
	auto res = DomainDiscrete_Div(alpha_1, M_1, t_div, alpha_2);
	alpha_1 = get<0>(res)[0];
	alpha_2 = get<1>(res)[0];
	double x_left = get<0>(res)[1];
	double x_right = get<0>(res)[2];
	int N_left = int(get<0>(res)[3]);
	int N = N_left + int(get<0>(res)[4]);
	double delta_x = get<0>(res)[5];
	double delta_tao_1 = get<0>(res)[6];
	double delta_tao_2 = get<1>(res)[6];
	int M_2 = int(get<1>(res)[7]);


	//change of variables:
	double a = r / pow(sigma, 2) - 0.5;
	double b = pow(a + 1, 2);
	double tao_div = 0.5*(T - t_div)*pow(sigma, 2);
	double tao_final = 0.5*T * pow(sigma, 2);

	//construct the result for each point
	mat Result_1(M_1 + 1, N + 1);
	Result_1.setZero();

	//calculte for tao = 0 (m=0)
	for (int i = 0; i < N + 1; i++) {
		double x_i = x_left + i*delta_x;
		Result_1(0, i) = K*exp(a*x_i)*max(exp(x_i) - 1, 0.0);
	}

	//calculate for x_left and x_right at [0, tao_div]
	for (int i = 0; i < M_1 + 1; i++) {
		double tao_i = i*delta_tao_1;
		//x_left
		Result_1(i, 0) = 0;
		//x_right
		Result_1(i, N) = K*exp(a*x_right + b*tao_i)*(exp(x_right) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}

	//construct maxtrix A(size: N-1)
	mat A_1(N - 1, N- 1);
	A_1.setZero();
	for (int i = 1; i < N- 2; i++) {
		A_1(i, i) = 1 + alpha_1;
		A_1(i, i - 1) = -0.5 * alpha_1;
		A_1(i, i + 1) = -0.5 * alpha_1;
	}
	A_1(0, 0) = 1 + alpha_1;
	A_1(0, 1) = -0.5 * alpha_1;
	A_1(N - 2, N - 2) = 1 + alpha_1;
	A_1(N - 2, N - 3) = -0.5 * alpha_1;

	//Get Tridiagonal_LU decomposition of A_1
	LSolver linearsolver;
	auto LU = linearsolver.tridiag_lu_no_pivoting(A_1);
	mat L = get<0>(LU);
	mat U = get<1>(LU);
	

	//calculte for m = 1:M_1 [0, tao_div]
	for (int m = 1; m < M_1 + 1; m++) {
		//construct vector b(size:N-1);
		vec b(N - 1);
		for (int i = 1; i < N - 2; i++) {
			b(i) = 0.5*alpha_1*Result_1(m - 1, i + 2) + (1 - alpha_1)*Result_1(m - 1, i + 1) + 0.5*alpha_1*Result_1(m - 1, i);
		}
		b(0) = 0.5*alpha_1*Result_1(m - 1, 2) + (1 - alpha_1)*Result_1(m - 1, 1)
			+ 0.5*alpha_1*Result_1(m - 1, 0) + 0.5*alpha_1*Result_1(m, 0);
		b(N - 2) = 0.5*alpha_1*Result_1(m - 1, N) + (1 - alpha_1)*Result_1(m - 1, N - 1)
			+ 0.5*alpha_1*Result_1(m - 1, N - 2) + 0.5*alpha_1*Result_1(m, N);

		//Solve for LUx=b : find y s.t. Ly=b (size of y: N-1)
		vec Y(N - 1);
		Y(0) = b(0);
		for (int i = 1; i < N - 1; i++) {
			Y(i) = b(i) - L(i, i - 1) * Y(i - 1);
		}

		//Solve for Ux=y £º find x s.t. Ux= y(size of x: N-1)
		Result_1(m, N - 1) = Y(N - 2) / U(N - 2, N - 2);
		for (int j = N - 3; j >= 0; j--) {
			Result_1(m, j + 1) = (Y(j) - U(j, j + 1) * Result_1(m, j + 2)) / U(j, j);
		}
	}

	//construct the result for each point
	mat Result_2(M_2 + 1, N + 1);
	Result_2.setZero();
	
	//calculate for x_left and x_right at [tao_div, tao_final]
	for (int i = 1; i < M_2 + 1; i++) {
		double tao_i = tao_div + i*delta_tao_2;
		//x_left
		Result_2(i, 0) = 0;
		//x_right
		Result_2(i, N) = K*exp(a*(x_right - log(1 - q)) + b*tao_i)*(exp(x_right) - exp(-2 * r*tao_i / pow(sigma, 2)));
	}

	//calculte for tao = tao_div (m=0)
	for (int i = 0; i < N + 1; i++) {
		Result_2(0, i) = pow(1 - q, -a)*Result_1(M_1, i);
	}


	//construct maxtrix A(size: N-1)
	mat A_2(N - 1, N - 1);
	A_2.setZero();
	for (int i = 1; i < N - 2; i++) {
		A_2(i, i) = 1 + alpha_2;
		A_2(i, i - 1) = -0.5 * alpha_2;
		A_2(i, i + 1) = -0.5 * alpha_2;
	}
	A_2(0, 0) = 1 + alpha_2;
	A_2(0, 1) = -0.5 * alpha_2;
	A_2(N - 2, N - 2) = 1 + alpha_2;
	A_2(N - 2, N - 3) = -0.5 * alpha_2;

	//Get Tridiagonal_LU decomposition of A_2
	LU = linearsolver.tridiag_lu_no_pivoting(A_2);
	L = get<0>(LU);
	U = get<1>(LU);

	//calculate for [tao_div, tao_final]
	for (int m = 1; m <= M_2; m++) {
		//construct vector b(size:N-1);
		vec b(N - 1);
		for (int i = 1; i < N - 2; i++) {
			b(i) = 0.5*alpha_2*Result_2(m - 1, i + 2) + (1 - alpha_2)*Result_2(m - 1, i + 1) + 0.5*alpha_2*Result_2(m - 1, i);
		}
		b(0) = 0.5*alpha_2*Result_2(m - 1, 2) + (1 - alpha_2)*Result_2(m - 1, 1)
			+ 0.5*alpha_2*Result_2(m - 1, 0) + 0.5*alpha_2*Result_2(m, 0);
		b(N - 2) = 0.5*alpha_2*Result_2(m - 1, N) + (1 - alpha_2)*Result_2(m - 1, N - 1)
			+ 0.5*alpha_2*Result_2(m - 1, N - 2) + 0.5*alpha_2*Result_2(m, N);

		//Solve for LUx=b : find y s.t. Ly=b (size of y: N-1)
		vec Y(N - 1);
		Y(0) = b(0);
		for (int i = 1; i < N - 1; i++) {
			Y(i) = b(i) - L(i, i - 1) * Y(i - 1);
		}

		//Solve for Ux=y £º find x s.t. Ux= y(size of x: N-1)
		Result_2(m, N - 1) = Y(N - 2) / U(N - 2, N - 2);
		for (int j = N - 3; j >= 0; j--) {
			Result_2(m, j + 1) = (Y(j) - U(j, j + 1) * Result_2(m, j + 2)) / U(j, j);
		}
	}
	
	
	//Pointwise convergence
	double x_compute = log(S / K);
	double x_1 = x_compute - delta_x;
	double x0 = x_compute;
	double x1 = x_compute + delta_x;

	double S_1 = K*exp(x_1);
	double S0 = K*exp(x0);
	double S1 = K*exp(x1);

	double V_1 = exp(-a*x_1 - b*tao_final)*Result_2(M_2, N_left - 1);
	double V0 = exp(-a*x0 - b*tao_final)*Result_2(M_2, N_left);
	double V1 = exp(-a*x1 - b*tao_final)*Result_2(M_2, N_left + 1);

	double Delta = (V1 - V_1) / (S1 - S_1);
	double Gamma = 2 * ((S0 - S_1)*V1 - (S1 - S_1)*V0 + (S1 - S0)*V_1) / ((S0 - S_1)*(S1 - S0)*(S1 - S_1));

	//Not sure to use delta_tao2 or not
	double delta_t = 2 * delta_tao_2 / pow(sigma, 2);
	double V_appro = exp(-a*x0 - b*(tao_final - delta_tao_2))*Result_2(M_2 - 1, N_left);
	double Option_Value = exp(-a*x_compute - b*tao_final)*Result_2(M_2, N_left);
	double Theta = (Option_Value - V_appro) / delta_t;


	vector<double> VG = { Result_2(M_2, N_left), Option_Value, Delta, Gamma, Theta };
	return make_tuple(VG, Result_1, Result_2);
}