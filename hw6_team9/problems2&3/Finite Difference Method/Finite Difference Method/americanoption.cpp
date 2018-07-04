#include <iostream>
#include <cmath>
#include "BlackScholes.h"
#include "Option.h"
#include "RWcsv.h"
#include <functional>
#include "BinomialTreePricer.h"
#include "equationsolver.h"
#include <tuple>

void AmericanPut()
{
	//Option(S0,K,T,sigma,q,r)
	Option op{ 42,40,0.75,0.32,0.02,0.04 };
	//Exact Value
	double exact = 3.3045362802172642;

	std::vector<double> bs = BSOption(op, 'p');
	Pricer amerput(op, 'p');

	//FD
	std::cout << "\nFinite Difference\n" << std::endl;
	int M = 4;
	Eigen::MatrixXd result(4, 9);
	std::vector<double> price;
	std::vector<double> error;
	std::vector<double> Euro;

	//double error = 0;
	for (int i = 0; i <4; i++)
	{
		price = amerput.AmericanOptionFiniteDifference(M*std::pow(4,i), 0.45);
		Euro.push_back(amerput.EuropreanOptionFiniteDifference(M*std::pow(4, i), 0.45));
		result(i, 0) = std::abs(price[0] - exact);
		result(i, 2) = std::abs(price[1] - exact);
		result(i, 4) = price[2];
		result(i, 5) = price[3];
		result(i, 6) = price[4];
		result(i, 7) = price[0] + (bs[0] - Euro[i]);
		result(i, 8) = std::abs(result(i, 7) - exact);
		price.clear();
	}
	for (int i = 1; i < 4; i++)
	{
		result(i, 1) = result(i, 0)/ result(i - 1, 0);
		result(i, 3) = result(i, 2)/ result(i - 1, 2);
	}

	save_csv(result, "american.csv");

}

void ImpliedVol()
{
	Option op{ 43,45,0.75,0,0.02,0.045 };
	BinomialTree bt(op);
	double market = 3.85;
	std::function<double(double)> binomialvol = [&](double vol) {op.sigma = vol; BinomialTree bt(op); return (bt.BTOptionPricer(2500, 'a', 'p'))[0] - market; };
	auto vol = Secant(0.1, 0.5, binomialvol);
	std::cout << std::get<0>(vol) << std::endl;
	std::cout << std::get<1>(vol) << std::endl;


}
int main()
{
	AmericanPut();
	ImpliedVol();
	return 0;
}