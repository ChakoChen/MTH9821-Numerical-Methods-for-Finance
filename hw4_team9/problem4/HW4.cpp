#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include "RWcsv.h"
#include "RandomNumberGenerator.h"
#include "BlackScholes.h"

using namespace Eigen;


int n = 10000*std::pow(2,9);
int seed = 1;
//(k,s0,q,vol,r,t)
Option op = { 55, 50, 0, 0.3, 0.04, 0.75 };
Option op1 = { 50,25,0,0.3,0.05,0.75 };
Option op2 = { 50,30,0,0.2,0.05,0.75 };

void LCtest()
{
	std::vector<double> unif=NLinearCongruential(seed, n);
	VectorXd result = VecTransform<double>(unif);
	save_csv<VectorXd>(result, "LCtest.csv");
}

void ITtest()
{
	std::vector<double> it = InverseTransform(seed, n);
	VectorXd result = VecTransform<double>(it);
	save_csv<VectorXd>(result, "ITtest.csv");
}
void ARtest()
{
	std::vector<double> ar = AcceptanceRejection(seed, n);
	VectorXd result = VecTransform<double>(ar);
	save_csv<VectorXd>(result, "ARtest.csv");
}

void BMtest()
{
	std::vector<double> bm = BoxMuller(seed, n);
	VectorXd result = VecTransform<double>(bm);
	save_csv<VectorXd>(result, "BMtest.csv");
}

double max(double a, double b)
{//return maximum of a and b
	if (a > b) return a;
	else
		return b;
}

double sum(std::vector<double> src)
{//return summation of entries of a vector
	double s = 0;
	for_each(src.begin(), src.end(), [&](double x) {s = s + x; });
	return s;
}

//Calculate Black-scholes value V_bs of the option
std::vector<double> BS = BSOption(op, 'p');

void ControlVariate()
{
	std::vector<double> z = BoxMuller(seed, n);

	std::vector<double> SpotPrice;	//S_i, i=0:n-1
	std::vector<double> Value;			//S_i, i=0:n-1
	double avgs, avgv;
	double bnom = 0, bdenom = 0;
	double b;
	std::vector<double> CVValue;
	std::vector<double> Vcv;

	int N=0;

	for (int k = 0; k < 10; k++)
	{
		N = 10000 * std::pow(2, k);

		std::cout << "N = " << N << " starts" << std::endl;

		SpotPrice.resize(N);
		Value.resize(N);
		CVValue.resize(N);
		bnom = 0;
		bdenom = 0;

		for (int i =0; i < N; i++)
		{
			SpotPrice[i] = op.S0*std::exp((op.r - op.vol*op.vol / 2.0)*op.T + op.vol*std::sqrt(op.T)*z[i]);
			Value[i] = std::exp(-op.r*op.T)*max(op.K - SpotPrice[i], 0);
		}

		avgs = sum(SpotPrice) / N;
		avgv = sum(Value) / N;

		for (int j=0; j < N; j++)
		{
			bnom = bnom + (SpotPrice[j] - avgs)*(Value[j] - avgv);
			bdenom = bdenom+std::pow((SpotPrice[j] - avgs), 2);
		}

		b = bnom / bdenom;

		for (int j = 0; j < N; j++)
		{
			CVValue[j] = Value[j] - b*(SpotPrice[j] - std::exp(op.r*op.T)*op.S0);
		}
		Vcv.push_back(sum(CVValue) / N);

	}

	MatrixXd result(10, 2);
	for (int j = 0; j < 10; j++)
	{
		result(j, 0) = Vcv[j];
		result(j, 1) = std::abs(Vcv[j] - BS[0]);
	}

	save_csv<MatrixXd>(result, "ControlVariate.csv");

	std::cout << "Finish!!" << std::endl;

}

void AntitheticVariables()
{
	std::vector<double> z1 = BoxMuller(seed, n/2);

	double spot1 = 0, spot2 = 0;
	double v1 = 0, v2 = 0;
	std::vector<double> SpotPrice;	//S_i, i=0:n-1
	std::vector<double> Value;			//S_i, i=0:n-1
	std::vector<double> Vav;

	int N = 0;

	for (int k = 0; k < 10; k++)
	{
		N = 10000 * std::pow(2, k);
		std::cout << "N = " << N << " starts" << std::endl;


		for (int i = 0; i < N/2; i++)
		{
			spot1 = op.S0*std::exp((op.r - op.vol*op.vol / 2.0)*op.T + op.vol*std::sqrt(op.T)*z1[i]);
			spot2 = op.S0*std::exp((op.r - op.vol*op.vol / 2.0)*op.T + op.vol*std::sqrt(op.T)*(-z1[i]));
			v1 = std::exp(-op.r*op.T)*max(op.K - spot1, 0);
			v2 = std::exp(-op.r*op.T)*max(op.K - spot2, 0);
			SpotPrice.push_back(spot1);
			SpotPrice.push_back(spot2);
			Value.push_back(v1);
			Value.push_back(v2);
		}
		Vav.push_back(sum(Value) / N);

		SpotPrice.clear();
		Value.clear();
	}

	MatrixXd result(10, 2);
	for (int j = 0; j < 10; j++)
	{
		result(j, 0) = Vav[j];
		result(j, 1) = std::abs(Vav[j] - BS[0]);
	}

	save_csv<MatrixXd>(result, "AntitheticVariables.csv");

	std::cout << "Finish!!" << std::endl;
}

void MomentMatching()
{
	std::vector<double> z = BoxMuller(seed, n);
	std::vector<double> SpotPrice;	//S_i, i=0:n-1
	std::vector<double> SpotPrice1;	//SpotPrice1 = SpotPrice *e^rT*S0/avgs
	std::vector<double> Value;			//S_i, i=0:n-1
	std::vector<double> Vmm;
	double avgs;
	int N = 0;

	for (int k = 0; k < 10; k++)
	{
		N = 10000 * std::pow(2, k);
		std::cout << "N = " << N << " starts" << std::endl;

		SpotPrice.resize(N);
		SpotPrice1.resize(N);
		Value.resize(N);

		for (int i = 0; i < N; i++)
		{
			SpotPrice[i] = op.S0*std::exp((op.r - op.vol*op.vol / 2.0)*op.T + op.vol*std::sqrt(op.T)*z[i]);
		}

		avgs = sum(SpotPrice) / N;

		for (int i = 0; i < N; i++)
		{
			SpotPrice1[i] = SpotPrice[i]*(std::exp(op.r*op.T)*op.S0)/avgs;
			Value[i] = std::exp(-op.r*op.T)*max(op.K - SpotPrice1[i], 0);
		}
		Vmm.push_back(sum(Value) / N);

		SpotPrice.clear();
		Value.clear();
	}

	MatrixXd result(10, 2);
	for (int j = 0; j < 10; j++)
	{
		result(j, 0) = Vmm[j];
		result(j, 1) = std::abs(Vmm[j] - BS[0]);
	}

	save_csv<MatrixXd>(result, "MomentMatching.csv");

	std::cout << "Finish!!" << std::endl;
}

void SimultaneousMMandCV()
{
	std::vector<double> z = BoxMuller(seed, n);

	std::vector<double> SpotPrice;	//S_i, i=0:n-1
	std::vector<double> SpotPrice1;	//SpotPrice1 = SpotPrice *e^rT*S0/avgs

	std::vector<double> Value;			//S_i, i=0:n-1
	std::vector<double> W;
	std::vector<double> Vcvmm;

	double avgs,avgv,avgw;
	double bnom = 0, bdenom = 0;
	double b;


	int N = 0;

	for (int k = 0; k < 10; k++)
	{
		N = 10000 * std::pow(2, k);
		std::cout << "N = " << N << " starts" << std::endl;

		SpotPrice.resize(N);
		SpotPrice1.resize(N);
		Value.resize(N);
		W.resize(N);

		for (int i = 0; i < N; i++)
		{
			SpotPrice[i] = op.S0*std::exp((op.r - op.vol*op.vol / 2.0)*op.T + op.vol*std::sqrt(op.T)*z[i]);	//S_i
		}

		avgs = sum(SpotPrice) / N;	//^S(n)

		for (int i = 0; i < N; i++)
		{
			SpotPrice1[i] = SpotPrice[i] * (std::exp(op.r*op.T)*op.S0) / avgs;
			Value[i] = std::exp(-op.r*op.T)*max(op.K - SpotPrice1[i], 0);
		}

		avgv = sum(Value) / N;

		for (int j = 0; j < N; j++)
		{
			bnom = bnom + (SpotPrice1[j] - std::exp(op.r*op.T)*op.S0)*(Value[j] - avgv);
			bdenom = bdenom + std::pow((SpotPrice[j] - std::exp(op.r*op.T)*op.S0), 2);
		}

		b = bnom / bdenom;

		for (int j = 0; j < N; j++)
		{
			W[j] = Value[j] - b*(SpotPrice1[j] - std::exp(op.r*op.T)*op.S0);
		}

		Vcvmm.push_back(sum(W) / N);

		SpotPrice.clear();
		SpotPrice1.clear();
		Value.clear();
		W.clear();
	}

	MatrixXd result(10, 2);
	for (int j = 0; j < 10; j++)
	{
		result(j, 0) = Vcvmm[j];
		result(j, 1) = std::abs(Vcvmm[j] - BS[0]);
	}

	save_csv<MatrixXd>(result, "SimultaneousMMandCV.csv");

	std::cout << "Finish!!" << std::endl;

}

void BasketOptions()
{
	//Generate 2n standard normal random variables
	int n1= 10000 * std::pow(2, 8);
	std::vector<double> z = BoxMuller(seed, 2*n1);
	std::vector<double> SpotPrice1;
	std::vector<double> SpotPrice2;
	std::vector<double> Value;
	std::vector<double> Vbo;
	double spot1, spot2, v;
	double corr = 0.25;

	int N = 0;
	for (int k = 0; k < 9; k++)
	{
		N = 10000 * std::pow(2, k);
		std::cout << "N = " << N << " starts" << std::endl;
		SpotPrice1.resize(N);
		SpotPrice2.resize(N);
		Value.resize(N);

		//calculate spot prices for stock1 and stock2
		for (int i = 0; i < N; i++)
		{
			spot1 = op1.S0*std::exp((op1.r - op1.vol*op1.vol / 2.0)*op1.T + op1.vol*std::sqrt(op1.T)*z[2 * i]);
			spot2 = op2.S0*std::exp((op2.r - op2.vol*op2.vol / 2.0)*op2.T + op2.vol*std::sqrt(op2.T)*(corr*z[2*i]+std::sqrt(1-corr*corr)*z[2*i+1]));
			SpotPrice1[i] = spot1;
			SpotPrice2[i] = spot2;
		}

		//Calculate the value of the basket option
		for (int i = 0; i < N; i++)
		{
			v = std::exp(-op1.r*op1.T)*max(SpotPrice1[i] + SpotPrice2[i]-op1.K, 0);
			Value[i] = v;
		}
		Vbo.push_back(sum(Value) / N);

		SpotPrice1.clear();
		SpotPrice2.clear();
		Value.clear();
	}

	MatrixXd result(9, 1);
	for (int j = 0; j < 9; j++)
	{
		result(j, 0) = Vbo[j];
	}

	save_csv<MatrixXd>(result, "BasketOptions.csv");

	std::cout << "Finish!!" << std::endl;

}

void PathBasketOption()
{
	//Generate 2n standard normal random variables
	int n1 = 50 * std::pow(2, 9);
	int m = 150;
	std::vector<double> z = BoxMuller(seed, 2 * n1*m);
	std::vector<double> Value;
	std::vector<double> Vpbo;
	double spot1, spot2, v;
	double corr = 0.25;
	double dt = op1.T / m;
	double maxpay = 0;
	int cnt=0;
	int N = 0;
	for (int k = 0; k < 10; k++)
	{
		N = 50 * std::pow(2, k);
		std::cout << "N = " << N << " starts" << std::endl;
		Value.resize(N);
		cnt = 0;

		for (int j = 0; j < N; j++)
		{
			spot1 = op1.S0;
			spot2 = op2.S0;
			maxpay = spot1 + spot2;

			for (int i = 0+m*cnt; i < m+m*cnt; i++)
			{
				spot1 = spot1*std::exp((op1.r - op1.vol*op1.vol / 2.0)*dt + op1.vol*std::sqrt(dt)*z[2 * i]);
				spot2 = spot2*std::exp((op2.r - op2.vol*op2.vol / 2.0)*dt + op2.vol*std::sqrt(dt)*(corr*z[2 * i] + std::sqrt(1 - corr*corr)*z[2 * i + 1]));
				maxpay = max(maxpay, spot1 + spot2);
			}
			cnt++;
			v = std::exp(-op1.r*op1.T)*max(maxpay-op1.K, 0);
			Value[j] = v;
		}
		Vpbo.push_back(sum(Value) / N);

		Value.clear();

	}

	MatrixXd result(10, 1);
	for (int j = 0; j < 10; j++)
	{
		result(j, 0) = Vpbo[j];
	}

	save_csv<MatrixXd>(result, "PathBasketOption.csv");

	std::cout << "Finish!!" << std::endl;

}

//Extra
void HestonModel()
{
	double s0 = 50;
	double vol = 0.3;
	double V0 = 0.09;
	double Var,mVar;
	double lambda = 3;
	double ltvar = std::pow(0.35,2);	//long term standard-deviation
	double sd = 0.25;	//standard deviation of the variance
	double corr = -0.15;
	double T = 0.5;
	double K = 50;
	double r = 0.05;

	//(k,s0,q,vol,r,t)
	Option op3 = { K,s0,0,vol,r,T };

	//time interval
	int m = 175;
	int n2 = 500 * std::pow(2, 5);

	std::vector<double> z = BoxMuller(seed, 2 * n2*m);
	std::vector<double> Value;
	std::vector<double> Vh;
	double approx = 0;
	std::vector<double> impvol;
	double spot, v;
	double dt = T / m;

	int cnt = 0;
	int N = 0;
	for (int k = 0; k < 6; k++)
	{
		N = 500 * std::pow(2, k);
		std::cout << "N = " << N << " starts" << std::endl;
		Value.resize(N);
		cnt = 0;

		for (int j = 0; j < N; j++)
		{
			spot = s0;
			Var=V0;		//variance	

			for (int i = 0 + m*cnt; i < m + m*cnt; i++)
			{
				mVar = max(Var, 0);
				Var = mVar - lambda*(mVar - ltvar)*dt + sd*std::sqrt(mVar)*std::sqrt(dt)*(corr*z[2 * i] + std::sqrt(1 - corr*corr)*z[2 * i + 1]);
				spot = spot*std::exp((r - mVar / 2.0)*dt + std::sqrt(mVar)*std::sqrt(dt)*z[2 * i]);
			}
			cnt++;
			v = std::exp(-r*T)*max(K-spot, 0);
			Value[j] = v;
		}
		approx = sum(Value) / N;
		Vh.push_back(approx);
		impvol.push_back(BSImpliedVol(op3, vol + 0.1, vol + 0.2, approx, 'p'));

		Value.clear();

	}

	MatrixXd result(6, 2);
	for (int j = 0; j < 6; j++)
	{
		result(j, 0) = Vh[j];
		result(j, 1) = impvol[j];
	}

	save_csv<MatrixXd>(result, "HestonModel.csv");

	std::cout << "Finish!!" << std::endl;


}

int main()
{
	//LCtest();
	//ITtest();
	//ARtest();
	//BMtest();

	ControlVariate();

	AntitheticVariables();
	
	MomentMatching();

	SimultaneousMMandCV();

	BasketOptions();

	PathBasketOption();

	HestonModel();
	return 0;
}
