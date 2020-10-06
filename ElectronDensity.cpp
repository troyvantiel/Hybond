#include "EDfunc.hpp"
#include <iostream>
#include <vector>
#include <math.h>

#define PI 3.14159265
#define abso(x) ((x > 0)? x: -x)			//defining the variables needed for the calculations
#define A(x,y) (elecHess[y * 3 + x - 4])

const double rhoCutoff = 0.1; //rho cut off variable
const int ArCutOff = 18; //argon cut off as the approximation only goes this far


using namespace std;

double power(double b, int p)
{
	double result =1;
	for(; p; --p)
	{
		result *=b;
	}
	return result;
}


void setup()
{
	czeta1 = new double[ArCutOff];
	czeta2 = new double[AtCutOff];
	czeta3 = new double[ArCutOff];  //creates the arrays with the argon cutoff size
	czzeta1 = new double[ArCutOff];
	czzeta2 = new double[ArCutOff];
	czzeta3 = new double[ArCutOff];
	elecHess = new double[9];
}
void getHesRho()
{
	rho =0;
	dxrho = 0; //set up variables with a rho and a rho for x,y and z
	dyrho = 0;
	dzrho = 0;
	for(int i = 0; i<9; i++)
	{
		elecHess[i] = 0; //fill up the elecHess array with 0's

	}
	for(int i =0; i <centers; i++)
	{
		int type = centertypes[i];
		double r = sqrt(distanceFromCenter[i]);
		double r1 = 1/r;
		double exp1 = exp(-r/zeta1[type]);
		double exp2 = exp(-r/zeta2[type]);
		double exp3 = exp(-r/zeta3[type]);
		double fac0 = c1[type] * exp1 + c2[type] * exp2 + c3[type] * exp3;
		double fac1 = czeta1[type] * exp1 + czeta2[type] * exp2 + czeta3[type] * exp3;
		double fac2 = czzeta1[type] * exp1 + czzeta2[type] * exp2 + czzeta3[type] * exp3;
		double fac3 = fac2 + fac1 * r1;


	}
}
double RDGrho(double x, double y, double z, double *srho)
{

}









