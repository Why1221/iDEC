#include <iostream>
#include <math.h>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "ChiSquared.h"

#define DBL_MAX 1e20
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

double chi2inv(double x, int n)
{
	return gaminv(x, n/2.0);
}

double gaminv(double x,double a)
{
	double y,y_new,y_old,h,cdf,pdf,h_abs;

	double eps=2.2204e-16;	

		if(x<eps)
		{
			y=sqrt(eps);
		}
		else
		{
			y=a*2;
		}

		y_old=y;

	int count;
	for(count=0;count<100;count++)
	{
		cdf=gamcdf(y_old,a);
		pdf=gampdf(y_old,a);
		
			h=(cdf-x)/pdf;
			y_new=y_old-h;
			if(y_new<=eps)
			{
				y_new=y_old/10;
				h=y_old-y_new;
			}

		h_abs=fabs(h);

		if(h_abs<sqrt(eps))goto loop;

			y_old=y_new;

	}
loop:

	return y_new;
}

double gampdf(double x,double a)
{
	double pdf;

	if(x>0)
	{
		pdf=exp(-a*log(2.0)+(a-1)*log(x)-x/2.0-LogGamma(a));	
	}
	else
	{
		pdf=0;
	}

	return pdf;
}

double gamcdf(double x,double a)
{
	double cdf;


	if(x!=0)
	{
		cdf=kf_gammap((double)a,x/2.0);
	}
	else
	{
		cdf=0;
	}



	return cdf;
}

double Gamma
	(
	double x    // We require x > 0
	)
{
	if (x <= 0.0)
	{
		std::stringstream os;
		os << "Invalid input argument " << x <<  ". Argument must be positive.";
		throw std::invalid_argument( os.str() ); 
	}

	// Split the function domain into three intervals:
	// (0, 0.001), [0.001, 12), and (12, infinity)

	///////////////////////////////////////////////////////////////////////////
	// First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

	if (x < 0.001)
		return 1.0/(x*(1.0 + gamma*x));

	///////////////////////////////////////////////////////////////////////////
	// Second interval: [0.001, 12)

	if (x < 12.0)
	{
		// The algorithm directly approximates gamma over (1,2) and uses
		// reduction identities to reduce other arguments to this interval.

		double y = x;
		int n = 0;
		bool arg_was_less_than_one = (y < 1.0);

		// Add or subtract integers as necessary to bring y into (1,2)
		// Will correct for this below
		if (arg_was_less_than_one)
		{
			y += 1.0;
		}
		else
		{
			n = static_cast<int> (floor(y)) - 1;  // will use n later
			y -= n;
		}

		// numerator coefficients for approximation over the interval (1,2)
		static const double p[] =
		{
			-1.71618513886549492533811E+0,
			2.47656508055759199108314E+1,
			-3.79804256470945635097577E+2,
			6.29331155312818442661052E+2,
			8.66966202790413211295064E+2,
			-3.14512729688483675254357E+4,
			-3.61444134186911729807069E+4,
			6.64561438202405440627855E+4
		};

		// denominator coefficients for approximation over the interval (1,2)
		static const double q[] =
		{
			-3.08402300119738975254353E+1,
			3.15350626979604161529144E+2,
			-1.01515636749021914166146E+3,
			-3.10777167157231109440444E+3,
			2.25381184209801510330112E+4,
			4.75584627752788110767815E+3,
			-1.34659959864969306392456E+5,
			-1.15132259675553483497211E+5
		};

		double num = 0.0;
		double den = 1.0;
		int i;

		double z = y - 1;
		for (i = 0; i < 8; i++)
		{
			num = (num + p[i])*z;
			den = den*z + q[i];
		}
		double result = num/den + 1.0;

		// Apply correction if argument was not initially in (1,2)
		if (arg_was_less_than_one)
		{
			// Use identity gamma(z) = gamma(z+1)/z
			// The variable "result" now holds gamma of the original y + 1
			// Thus we use y-1 to get back the orginal y.
			result /= (y-1.0);
		}
		else
		{
			// Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
			for (i = 0; i < n; i++)
				result *= y++;
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	// Third interval: [12, infinity)

	if (x > 171.624)
	{
		// Correct answer too large to display. Force +infinity.
		double temp = DBL_MAX;
		return temp*2.0;
	}

	return exp(LogGamma(x));
}

double LogGamma
	(
	double x    // x must be positive
	)
{
	if (x <= 0.0)
	{
		std::stringstream os;
		os << "Invalid input argument " << x <<  ". Argument must be positive.";
		throw std::invalid_argument( os.str() ); 
	}

	if (x < 12.0)
	{
		return log(fabs(Gamma(x)));
	}

	// Abramowitz and Stegun 6.1.41
	// Asymptotic series should be good to at least 11 or 12 figures
	// For error analysis, see Whittiker and Watson
	// A Course in Modern Analysis (1927), page 252

	static const double c[8] =
	{
		1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
	};
	double z = 1.0/(x*x);
	double sum = c[7];
	for (int i=6; i >= 0; i--)
	{
		sum *= z;
		sum += c[i];
	}
	double series = sum/x;

	static const double halfLogTwoPi = 0.91893853320467274178032973640562;
	double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return logGamma;
}

double kf_lgamma(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

static double _kf_gammap(double s, double z)
{
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}

static double _kf_gammaq(double s, double z)
{
	int j;
	double C, D, f;
	f = 1. + z - s; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for (j = 1; j < 100; ++j) {
		double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
		D = b + a * D;
		if (D < KF_TINY) D = KF_TINY;
		C = b + a / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s) - log(f));
}

double kf_gammap(double s, double z)
{
	return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);
}

double kf_gammaq(double s, double z)
{
	return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);
}
