#include <algorithm>
#include "BinomialPricing.hpp"
using namespace std;

vector<double> Binomial::European_Nodiv(string type, int N){
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp((r - q)*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	vector<double> Value;
	vector<double> Delta;
	vector<double> Gamma;

	//pay off of an European PUT option
	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value.push_back(max(0.0, K - S*pow(u, N - i)*pow(d, i)));
		}
	}
	//pay off of an European CALL option
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value.push_back(max(0.0, S*pow(u, N - i)*pow(d, i) - K));
		}
	}

	//calculate values in the middle
	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			Value[i] = p1*Value[i] + p2*Value[i + 1];
		}
		if (k == 1) {
			//record the value to calculate the Delta
			Delta = { Value[0], Value[1] };
		}
		else if (k == 2) {
			//record the value to calculate the Gamma
			Gamma = { Value[0],Value[1],Value[2] };
		}
	}
	
	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::American_NoDiv(string type, int N) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();
	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp((r - q)*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	vector<double> Value;
	vector<double> Delta;
	vector<double> Gamma;
	//pay off of an European PUT option
	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value.push_back(max(0.0, K - S*pow(u, N - i)*pow(d, i)));
		}
	}
	//pay off of an European CALL option
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value.push_back(max(0.0, S*pow(u, N - i)*pow(d, i) - K));
		}
	}

	//calculate value in the middle
	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			double earlyexercise;
			//American Put
			if (type.compare("PUT") == 0) {
				earlyexercise = max(0.0, K - S*pow(u, k - i)*pow(d, i));
			}
			//American Call
			else if (type.compare("CALL") == 0) {
				earlyexercise = max(0.0, S*pow(u, k - i)*pow(d, i) - K);
			}
			//compare the value of holding the option and early excersice
			double holding = p1*Value[i] + p2*Value[i + 1];
			Value[i] = max(holding, earlyexercise);
		}
		if (k == 1) {
			//Record the value to calculate the Delta
			Delta = { Value[0], Value[1] };
		}
		else if (k == 2) {
			//Record the value to calculate the Gamma
			Gamma = { Value[0],Value[1],Value[2] };
		}
	}
	
	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	//Gamma
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::European_PropDiv(string type, int N) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	vector<double> Value(N + 1, 0), Delta(2, 0), Gamma(3, 0);

	//pay off of an European PUT option
	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, K - (1 - q)*S*pow(u, N - i)*pow(d, i));
		}
	}
	//pay off of an European CALL option
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, (1 - q)*S*pow(u, N - i)*pow(d, i) - K);
		}
	}
	
	//calculate value in the middle
	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			//calculate the present value
			Value[i] = p1*Value[i] + p2*Value[i + 1];
		}
		if (k == 1) {
			//Record the value to calculate the Delta
			Delta = { Value[0], Value[1] };
		}
		else if (k == 2) {
			//Record the value to calculate the Gamma
			Gamma = { Value[0],Value[1],Value[2] };
		}
	}
	
	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::American_PropDiv(string type, int N, double t_div) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();
	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);


	int DIV = int(t_div / deltaT);
	vector<double> Value(N + 1, 0), Delta, Gamma;

	//American Put
	if (type.compare("PUT") == 0) {
		//calculate payoff of the leave nodes
		for (int i = 0; i <= N; i++) {
			Value[i] = max(K - (1 - q)*S*pow(u, N - i)*pow(d, i), 0.0);
		}

		//calculate the values in the middle
		for (int n = N - 1; n >= DIV; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], K - (1 - q)*S*pow(u, n - i)*pow(d, i));
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		for (int n = (DIV - 1); n >= 0; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], K - S*pow(u, n - i)*pow(d, i));
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
	}

	//American Call
	else if (type.compare("CALL") == 0) {
		//calculate payoff of the leave nodes
		for (int i = 0; i <= N; i++) {
			Value[i] = max((1 - q)*S*pow(u, N - i)*pow(d, i) - K, 0.0);
		}

		//calculate the values in the middle
		for (int n = N - 1; n >= DIV+1; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], (1 - q)*S*pow(u, n - i)*pow(d, i) - K);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		for (int n = DIV; n >= 0; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], S*pow(u, n - i)*pow(d, i) - K);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
	}

	vector<double> Result = { Value[0] };
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::European_PropDiv_Mult(string type, int N, int count) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	vector<double> Value(N + 1, 0), Delta(2, 0), Gamma(3, 0);
	double S_new = S*pow(1 - q, count);

	//pay off of an European PUT option
	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, K - S_new*pow(u, N - i)*pow(d, i));
		}
	}
	//pay off of an European CALL option
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, S_new*pow(u, N - i)*pow(d, i) - K);
		}
	}

	//calculate value in the middle
	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			//calculate the present value
			Value[i] = p1*Value[i] + p2*Value[i + 1];
		}
		if (k == 1) {
			//Record the value to calculate the Delta
			Delta = { Value[0], Value[1] };
		}
		else if (k == 2) {
			//Record the value to calculate the Gamma
			Gamma = { Value[0],Value[1],Value[2] };
		}
	}

	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;

}

vector<double> Binomial::American_PropDiv_Mult(string type, int N, vector<double> times) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);


	double t1 = times[0], t2 = times[1], t3 = times[2];
	int T1 = int(t1 / deltaT) + 1, T2 = int(t2 / deltaT) + 1, T3 = int(t3 / deltaT) + 1;
	double S_new = S*pow(1 - q, 3);

	vector<double> Value(N + 1, 0), Delta(2, 0), Gamma(3, 0);
	double early_exercise;
	if (type.compare("PUT") == 0) {
		//calculate payoff in the nodes
		for (int i = 0; i < N; i++) {
			Value[i] = max(K - S_new*pow(u, N - i)*pow(d, i), 0.0);
		}
		//calculate the values in the middle
		//[t3, T]
		for (int n = N - 1; n >= T3; n--) {
			S_new = S*pow(1 - q, 3);
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - S_new*pow(u, n - i)*pow(d, i));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
		}
		//[t2, t3]
		for (int n = T3 - 1; n >= T2; n--) {
			S_new = S*pow(1 - q, 2);
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - S_new*pow(u, n - i)*pow(d, i));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
		}
		//[t1, t2]
		for (int n = T2 - 1; n >= T1; n--) {
			S_new = S*pow(1 - q, 1);
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - S_new*pow(u, n - i)*pow(d, i));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
			if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		//[0, t1]
		for (int n = T1 - 1; n >= 0; n--) {
			S_new = S;
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - S_new*pow(u, n - i)*pow(d, i));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}

	}

	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::European_FixDiv(string type, int N, double t_div, double v_div) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	vector<double> Value(N+1, 0), Delta, Gamma;

	//pay off of an European PUT option
	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, K - (S - v_div*exp(-r*t_div))*pow(u, N - i)*pow(d, i));
		}
	}
	//pay off of an European CALL option
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, (S - v_div*exp(-r*t_div))*pow(u, N - i)*pow(d, i) - K);
		}
	}

	//calculate value in the middle
	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			//calculate the present value
			Value[i] = p1*Value[i] + p2*Value[i + 1];
		}
		if (k == 1) {
			//Record the value to calculate the Delta
			Delta = { Value[0], Value[1] };
		}
		else if (k == 2) {
			//Record the value to calculate the Gamma
			Gamma = { Value[0],Value[1],Value[2] };
		}
	}

	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::American_FixDiv(string type, int N, double t_div, double v_div) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);


	int DIV = int(t_div / deltaT);
	vector<double> Value(N + 1, 0), Delta(2), Gamma(3);

	//American Put
	if (type.compare("PUT") == 0) {
		//calculate payoff of the leave nodes
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, K - (S - v_div*exp(-r*t_div))*pow(u, N - i)*pow(d, i));
		}

		//calculate the values in the middle
		for (int n = N - 1; n >= DIV; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], K - (S - v_div*exp(-r*t_div))*pow(u, n - i)*pow(d, i));
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		for (int n = (DIV - 1); n >= 0; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], K - S*pow(u, n - i)*pow(d, i));
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
	}

	//American Call
	else if (type.compare("CALL") == 0) {
		//calculate payoff of the leave nodes
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, (S - v_div*exp(-r*t_div))*pow(u, N - i)*pow(d, i) - K);
		}

		//calculate the values in the middle
		for (int n = N - 1; n >= DIV + 1; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], (S - v_div*exp(-r*t_div))*pow(u, n - i)*pow(d, i) - K);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		for (int n = DIV; n >= 0; n--) {
			for (int i = 0; i <= n; i++) {
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], S*pow(u, n - i)*pow(d, i) - K);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
	}

	vector<double> Result = { Value[0] };
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::European_FixDiv_Mult(string type, int N, vector<double> times, double v_div) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();
	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	double t1 = times[0],t2 = times[1], t3 = times[2];

	double S_new = S - v_div*exp(-r*t1) - v_div*exp(-r*t2) - v_div*exp(-r*t3);

	vector<double> Value(N + 1, 0), Delta, Gamma;

	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, K - S_new*pow(u, N - i)*pow(d, i));
		}
	}
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, S_new*pow(u, N - i)*pow(d, i) - K);
		}
	}

	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			Value[i] = p1*Value[i] + p2*Value[i + 1];
		}
		if (k == 1) {
			Delta = { Value[0], Value[1] };
		}
		if (k == 2) {
			Gamma = { Value[0], Value[1], Value[2] };
		}
	}
	vector<double> Result = { Value[0] };
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::American_FixDiv_Mult(string type, int N, vector<double> times, double v_div) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);


	double t1 = times[0], t2 = times[1], t3 = times[2];
	int T1 = int(t1 / deltaT) + 1, T2 = int(t2 / deltaT) + 1, T3 = int(t3 / deltaT) + 1;

	double S_new = S - v_div*exp(-r*t1) - v_div*exp(-r*t2) - v_div*exp(-r*t3);

	vector<double> Value(N + 1, 0), Delta(2, 0), Gamma(3, 0);
	double early_exercise;
	if (type.compare("PUT") == 0) {
		//calculate payoff in the nodes
		for (int i = 0; i < N; i++) {
			Value[i] = max(K - S_new*pow(u, N - i)*pow(d, i), 0.0);
		}
		//calculate the values in the middle
		//[t3, T]
		for (int n = N - 1; n >= T3; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - S_new*pow(u, n - i)*pow(d, i));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
		}
		//[t2, t3]
		for (int n = T3 - 1; n >= T2; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - (S_new*pow(u, n - i)*pow(d, i) + v_div*exp(-r*t3 + r*deltaT*n)));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
		}
		//[t1, t2]
		for (int n = T2 - 1; n >= T1; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - (S_new*pow(u, n - i)*pow(d, i) + v_div*exp(-r*t3 + r*deltaT*n) + v_div*exp(-r*t2 + r*deltaT*n)));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
			if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		//[0, t1]
		for (int n = T1 - 1; n >= 0; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - (S_new*pow(u, n - i)*pow(d, i) + v_div*exp(-r*t3 + r*deltaT*n) + v_div*exp(-r*t2 + r*deltaT*n)
					+ v_div*exp(-r*t1 + r*deltaT*n)));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}

	}

	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::European_Complex(string type, int N, vector<double> dtime, vector<double> dvalue) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();
	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	double t1 = dtime[0];
	double t2 = dtime[1];
	double t3 = dtime[2];
	double v1 = dvalue[0];
	double v2 = dvalue[1];
	double v3 = dvalue[2];

	double S_new = (1 - v2)*(S - v1*exp(-r*t1)) - v3*exp(-r*t3);

	vector<double> Value(N + 1, 0), Delta, Gamma;

	if (type.compare("PUT") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, K - S_new*pow(u, N - i)*pow(d, i));
		}
	}
	else if (type.compare("CALL") == 0) {
		for (int i = 0; i <= N; i++) {
			Value[i] = max(0.0, S_new*pow(u, N - i)*pow(d, i) - K);
		}
	}

	for (int k = N - 1; k >= 0; k--) {
		for (int i = 0; i <= k; i++) {
			Value[i] = p1*Value[i] + p2*Value[i + 1];
		}
		if (k == 1) {
			Delta = { Value[0], Value[1] };
		}
		if (k == 2) {
			Gamma = { Value[0], Value[1], Value[2] };
		}
	}
	vector<double> Result = { Value[0] };
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}

vector<double> Binomial::American_Complex(string type, int N, vector<double> dtime, vector<double> dvalue) {
	double S = getS();
	double K = getK();
	double T = getT();
	double r = getr();
	double sigma = getsigma();
	double q = getq();

	double deltaT = T / N;
	double u = exp(sigma*sqrt(deltaT));
	double d = 1 / u;
	double p = (exp(r*deltaT) - d) / (u - d);		//risk-neutral prob. of price going up
	double p1 = exp(-r*deltaT)*p;
	double p2 = exp(-r*deltaT)*(1 - p);

	
	double t1 = dtime[0], t2 = dtime[1], t3 = dtime[2];
	int T1 = int(t1 / deltaT) + 1, T2 = int(t2 / deltaT) + 1, T3 = int(t3 / deltaT) + 1;

	double v1 = dvalue[0], v2 = dvalue[1], v3 = dvalue[2];

	double S_new = (1 - v2)*(S - v1*exp(-r*t1)) - v3*exp(-r*t3);

	vector<double> Value(N + 1, 0), Delta(2, 0), Gamma(3, 0);
	double early_exercise;
	if (type.compare("PUT") == 0) {
		//calculate payoff in the nodes
		for (int i = 0; i < N; i++) {
			Value[i] = max(K - S_new*pow(u, N - i)*pow(d, i), 0.0);
		}
		//calculate the values in the middle
		//[t3, T]
		for (int n = N - 1; n >= T3; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - S_new*pow(u, n - i)*pow(d, i));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
		}
		//[t2, t3]
		for (int n = T3 - 1; n >= T2; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - (S_new*pow(u, n - i)*pow(d, i) + v3*exp(-r*t3 + r*deltaT*n)));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
		}
		//[t1, t2]
		for (int n = T2 - 1; n >= T1; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - (S_new*pow(u, n - i)*pow(d, i) + v3*exp(-r*t3 + r*deltaT*n)) / (1 - v2));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
			if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
		//[0, t1]
		for (int n = T1 - 1; n >= 0; n--) {
			for (int i = 0; i <= n; i++) {
				early_exercise = max(0.0, K - ((S_new*pow(u, n - i)*pow(d, i) + v3*exp(-r*t3 + r*deltaT*n))/ (1 - v2)
					+ v1*exp(-r*t1 + r*deltaT*n)));
				Value[i] = max(p1*Value[i] + p2*Value[i + 1], early_exercise);
			}
			if (n == 1) {
				//Record the value to calculate the Delta
				Delta = { Value[0], Value[1] };
			}
			else if (n == 2) {
				//Record the value to calculate the Gamma
				Gamma = { Value[0],Value[1],Value[2] };
			}
		}
	}

	//Value Delta Gamma Theta 
	vector<double> Result;
	//Value 
	Result.push_back(Value[0]);
	//Delta
	Result.push_back((Delta[0] - Delta[1]) / (S*u - S*d));
	//Gamma
	double tempD1 = (Gamma[0] - Gamma[1]) / (S*u*u - S*u*d);
	double tempD2 = (Gamma[1] - Gamma[2]) / (S*u*d - S*d*d);
	Result.push_back(2 * (tempD1 - tempD2) / (S*u*u - S*d*d));
	//Theta 
	Result.push_back((Gamma[1] - Value[0]) / (2 * deltaT));

	return Result;
}
