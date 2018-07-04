#include <string>

class Option {
	//private menmbers
private:
	//interest rate, time horizon, strike price
	double r, T, K;
	//initial strock prices, volatility, divident
	double S0, sigma, q;
public:
	Option();
	//copy const
	Option(const Option &);
	//initializing
	Option(double S0, double T, double K, double r, double sigma, double q);
	virtual ~Option();
	//copy assign operator
	Option& operator = (const Option &);
	//get initial value S0;
	double SpotPrice();
	//get initial value of volatility
	double Volatility();
	//get time horizon
	double Maturity();
	//get strike price
	double Strike();
	//get dividend
	double Dividend();
	//get interest rate
	double Interest();
	//get option type
	virtual std::string OptionType() = 0;
	//get option exact value
	virtual double Price() = 0;
};
