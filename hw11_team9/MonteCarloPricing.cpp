#include "MonteCarloPricing.hpp"
#include "NumberGenerators.hpp"

double MeanofVector(vector<double> Input) {
	double sum = 0.0;
	int n = Input.size();
	for (int i = 0; i < n; i++) {
		sum += Input[i];
	}
	return (sum / n);
}

vector<double> MonteCarlo::PathDependent_Div(int N, vector<double> Div_times, vector<double> Div_values) {
	NumberGenerators NG;
	vector<double> Random = NG.BoxMuller(N);
	double S = getS();
	double T = getT();
	double K = getK();
	double r = getr();
	double sigma = getsigma();

	
	int n = N / 4;
	
	double A_Coeff = r - 0.5*sigma*sigma;
	double t1 = Div_times[0];
	double t2 = Div_times[1] - Div_times[0];
	double t3 = Div_times[2] - Div_times[1];
	double t4 = T - Div_times[2];

	double S_tlt, S_temp;
	vector<double> P, P_tlt, Delta, Delta_tlt;
	for (int i = 0; i < n; i++) {
		S_temp = S*exp(A_Coeff*t1 + sigma*Random[4 * i] * sqrt(t1)) - Div_values[0];
		S_temp = S_temp * exp(A_Coeff*t2 + sigma*Random[4 * i + 1] * sqrt(t2))*(1 - Div_values[1]);
		S_temp = S_temp * exp(A_Coeff*t3 + sigma*Random[4 * i + 2] * sqrt(t3)) - Div_values[2];
		S_temp = S_temp * exp(A_Coeff*t4 + sigma*Random[4 * i + 3] * sqrt(t4));
		
		S_tlt = S*exp(A_Coeff*T + sigma*(sqrt(t1)*Random[4 * i] + sqrt(t2)*Random[4 * i + 1]
			+ sqrt(t3)*Random[4 * i + 2] + sqrt(t4)*Random[4 * i + 3]));

		//Value
		P.push_back(exp(-r*T)*max(K - S_temp, 0.0));
		P_tlt.push_back(exp(-r*T)*max(K - S_tlt, 0.0));
		
		//Dalta
		Delta.push_back(K > S_temp ? -1 * exp(-r*T)*S_tlt / S* (1 - Div_values[1]) : 0);
		Delta_tlt.push_back(K > S_tlt ? -1 * exp(-r*T)*S_tlt / S : 0);
	}
	
	double P_mean = MeanofVector(P);
	double P_tlt_mean = MeanofVector(P_tlt);
	
	//b_value
	double sum1 = 0, sum2 = 0;
	for (int i = 0; i < n; i++) {
		sum1 += (P[i] - P_mean)*(P_tlt[i] - P_tlt_mean);
		sum2 += pow(P_tlt[i] - P_tlt_mean, 2);
	}
	double b = sum1 / sum2;
	
	//P_CV
	double P_CV = 0;
	OptionPricing Opt(S, K, T, r, 0, sigma);
	double P_BS = Opt.BlackScholes("PUT")[0];
	for (int i = 0; i < n; i++) {
		P_CV += P[i] - b*(P_tlt[i] - P_BS);
	}
	P_CV = P_CV / n;

	//Delta_Mean, Delta_tlt_Mean
	double Delta_Mean = MeanofVector(Delta);
	double Delta_tlt_Mean = MeanofVector(Delta_tlt);

	//b_delta
	sum1 = 0, sum2 = 0;
	for (int i = 0; i < n; i++) {
		sum1 += (Delta[i] - Delta_Mean)*(Delta_tlt[i] - Delta_tlt_Mean);
		sum2 += pow(Delta_tlt[i] - Delta_tlt_Mean, 2);
	}
	double b_delta = sum1 / sum2;

	//Delta_CV
	double Delta_CV = 0;
	double Delta_BS = Opt.BlackScholes("PUT")[1];
	for (int i = 0; i < n; i++) {
		Delta_CV += Delta[i] - b_delta*(Delta_tlt[i] - Delta_BS);
	}
	Delta_CV = Delta_CV / n;

	
	vector<double> res = { P_mean, Delta_Mean, P_CV, Delta_CV};

	return res;
}

