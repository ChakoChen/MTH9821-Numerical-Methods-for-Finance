#include "Option.hpp"

class AmericanPut : public Option
{
private:
	double V_exact;
    
public:
	AmericanPut();
	AmericanPut(const AmericanPut &);
	AmericanPut(double S0, double T, double K, double r, double sigma, double q);
	AmericanPut(double S0, double T, double K, double r, double sigma, double q, double exact);
	virtual ~AmericanPut();
	AmericanPut& operator = (const AmericanPut &);
	virtual std::string OptionType();
	virtual double Price();
};
