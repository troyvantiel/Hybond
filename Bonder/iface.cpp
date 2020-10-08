#include "stdafx.h"
#include "iface.h"
#include "fill.h"
#include <math.h>
#include <vector>

#define PI           3.14159265358979323846
#define abso(x) ((x > 0)? x: -x)
#define A(x,y) (elecHess[y * 3 +x - 4])

using namespace std;
const double rhoCutoff = 0.1;

//must be set up


int type2ix[56] = { 0, 1, 0, 0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5 };
int type2iy[56] = { 0, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 };
int type2iz[56] = { 0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 1, 0, 2, 2, 1, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0, 5, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0 };


//exponental function
inline double power(double b, int p)
{
	double result = 1;
	for (; p; --p)
	{
		result *= b;
	}
	return result;
}

//setuo
void analysisBatch::setUpBatch(double x, double y, double z, double resalution)
{
	cout << "setting up the batch" << x<< y<< z << endl;
	offsetx = x;
	offsety = y;
	offsetz = z;
	res = resalution;
}

//more setup
analysisBatch::analysisBatch(wfnData input)
{

	elecHess = new double[9];

	centers = input.nuc;
	

	distanceFromCenter = new double[centers];
	dx = new double[centers];
	dy = new double[centers];
	dz = new double[centers];
	
	centertypes = input.type;
	centerXvaule = input.x;
	centerYvaule = input.y;
	centerZvaule = input.z;
	czeta1 = new double[maxAtom];
	czeta2 = new double[maxAtom];
	czeta3 = new double[maxAtom];
	czzeta1 = new double[maxAtom];
	czzeta2 = new double[maxAtom];
	czzeta3 = new double[maxAtom];


	for (int i= 0; i < maxAtom; i++)
	{
		czeta1[i] = c1[i] / zeta1[i];
		czeta2[i] = c2[i] / zeta2[i];
		czeta3[i] = c3[i] / zeta3[i];
		czzeta1[i] = czeta1[i] / zeta1[i];
		czzeta2[i] = czeta2[i] / zeta2[i];
		czzeta3[i] = czeta3[i] / zeta3[i];
	}



}

//cleanup
analysisBatch::~analysisBatch()
{
	delete [] elecHess;
	delete [] distanceFromCenter;
	delete [] dx;
	delete [] dy;
	delete [] dz;
	delete [] czeta1;
	delete [] czeta2;
	delete [] czeta3;
	delete [] czzeta1;
	delete [] czzeta2;
	delete [] czzeta3;
	if (isCulled)
	{
		delete [] centerXvaule;
		delete [] centerYvaule;
		delete [] centerZvaule;
		delete [] centertypes;
	}
}


//clears MO wavefunction data
void analysisBatch::vectorReset(double x,double y, double z)
{
	for (size_t i = 0; i < centers; ++i)
	{
		dx[i] = (x - centerXvaule[i]);
		dy[i] = (y - centerYvaule[i]);
		dz[i] = (z - centerZvaule[i]);

		distanceFromCenter[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	}
}





double analysisBatch::atomx(int i)
{
	return centerXvaule[i];
}

double analysisBatch::atomy(int i)
{
	return centerYvaule[i];
}

double analysisBatch::atomz(int i)
{
	return centerZvaule[i];
}






double analysisBatch::getSignOfSecondEiganVaule()
{
	//get eiganVaules
	double eigan[3],phi;

	double p1 = A(1, 2) * A(1, 2) + A(1, 3) * A(1, 3) + A(2, 3) * A(2, 3);
	if (p1 == 0)
	{
		// A is diagonal.
		printf("does this happen?");
		eigan[0] = A(1, 1);
		eigan[1] = A(2, 2);
		eigan[2] = A(3, 3);
	}
	else
	{
		double q = (A(1,1) + A(2,2) + A(3,3)) / 3;
		double p2 = (A(1, 1) - q) * (A(1, 1) - q) + (A(2, 2) - q) * (A(2, 2) - q) + (A(3, 3) - q) * (A(3, 3) - q) + 2 * p1;
		double p = sqrt(p2 / 6);
		double B[9];
		for (size_t i = 0; i < 9; ++i)
		{
			B[i] = (1 / p) * (elecHess[i] - q * ((i % 4) ? 0 : 1)) ; // I is the identity matrix;
		}
		
		double r = (B[0] * B[4] * B[8] - B[0] * B[5] * B[7] - B[1] * B[3] * B[8] + B[1] * B[5] * B[6] + B[2] * B[3] * B[7] - B[2] * B[4] * B[6]) / 2;

		// In exact arithmetic for a symmetric matrix - 1 <= r <= 1;
		// but computation error can leave it slightly outside this range.
		if (r <= -1)
			phi = PI / 3;
		else if (r >= 1)
			phi = 0;
		else
			phi = acos(r) / 3;

		//the eigenvalues satisfy eig3 <= eig2 <= eig1
		eigan[0] = q + 2 * p * cos(phi);
		eigan[2] = q + 2 * p * cos(phi + (2 * PI / 3));
		eigan[1] = 3 * q - eigan[0] - eigan[2];    // since trace(A) = eig1 + eig2 + eig3;
	}
	//printf("%f\n",eigan[1]);
	return ((eigan[1] > 0) ? 1 : -1);
}

void analysisBatch::getDRho()
{
	//cout << "Starting calculation of DRho" << endl;
	rho = 0;
	dxrho = 0;
	dyrho=0;
	dzrho=0;
	for(int i = 0; i <centers; i++)
	{
		int type = centertypes[i];
		double r = sqrt(distanceFromCenter[i]);
		double r1 = 1/r;
		double exp1 = exp(-r/zeta1[type]);
		double exp2 = exp(-r/zeta2[type]);
		double exp3 = exp(-r/zeta3[type]);
		double fac0 = c1[type] * exp1 + c2[type] * exp2 + c3[type] * exp3;
		double fac1 = czeta1[type]*exp1 + czeta2[type]*exp2 + czeta3[type]*exp3;
		rho += fac0;
		dxrho -= fac1 * dx[i] * r1;
		dyrho -= fac1 * dy[i] * r1;
		dzrho -= fac1 * dz[i] * r1;
	}
}


void analysisBatch::getHesRho()
{
	//cout << "Starting calculation of HesRho" << endl;
	rho = 0;
	dxrho = 0;
	dyrho=0;
	dzrho=0;
	for (int i = 0 ; i < 9; i++)
		elecHess[i] = 0;

	for(int i = 0; i < centers; i++)
	{
		int type = centertypes[i];
		double r = sqrt(distanceFromCenter[i]);
		double r1 = 1/r;
		double exp1 = exp(-r/zeta1[type]);
		double exp2 = exp(-r/zeta2[type]);
		double exp3 = exp(-r/zeta3[type]);
		double fac0 = c1[type] * exp1 + c2[type] * exp2 + c3[type] * exp3;
		double fac1 = czeta1[type]*exp1 + czeta2[type]*exp2 + czeta3[type]*exp3;
		double fac2 = czzeta1[type]*exp1 + czzeta2[type]*exp2 + czzeta3[type]*exp3;
		double fac3 = fac2 + fac1 * r1;
		rho += fac0;
		dxrho -= fac1 * dx[i] * r1;
		dyrho -= fac1 * dy[i] * r1;
		dzrho -= fac1 * dz[i] * r1;
		elecHess[0] += fac3 * dx[i] * r1* dx[i] * r1 - fac1 * r1;
		elecHess[4] += fac3 * dy[i] * r1* dy[i] * r1 - fac1 * r1;
		elecHess[8] += fac3 * dz[i] * r1* dz[i] * r1 - fac1 * r1;
		elecHess[1] += r1 * r1 * dy[i] * dx[i] * fac3;
		elecHess[2] += r1 * r1 * dz[i] * dx[i] * fac3;
		elecHess[5] += r1 * r1 * dz[i] * dy[i] * fac3;

	}
	elecHess[3] = elecHess[1];
	elecHess[6] = elecHess[2];
	elecHess[7] = elecHess[5];
}


double analysisBatch::RDG_rho(double x, double y, double z,double *srho)
{

	vectorReset(x, y, z);
	getHesRho();
	double rhosq = dxrho*dxrho + dyrho*dyrho + dzrho*dzrho;

	double absrho = (rho < 0) ? -rho : rho;
	*srho = absrho * getSignOfSecondEiganVaule();

	if (rho == 0)
	{
		return 999;
	}




	if (absrho > rhoCutoff)
		return 999;

	return 0.161620459673995 * sqrt(rhosq) / (pow(rho, 4.0 / 3));
	//0.161620459673995D0 =  1/(2*(3*pi^2)^(1/3))
}

/*
   double analysisBatch::LOL(double x, double y, double z)
   {
   double rho = 0;
   double kineng = 0;
   for (size_t i = 0; i < nmo; ++i)
   {
   rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
   kineng += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionDZ[i] * moWavefuntionDZ[i]);
   }
   kineng /= 2;
   double Fc = 2.871234000; //magic number do not touch
   double Dh = Fc*pow(rho,(5.0 / 3.0));
   kineng = Dh / kineng;
   return 1.0/(1.0 + kineng);
   }
   */

double analysisBatch::elf()
{
	double Fc = 2.871234000;
	double rhosq = dxrho * dxrho + dyrho * dyrho + dzrho * dzrho;
	double lap = elecHess[0] + elecHess[4] + elecHess[8];
	double kineng = (3.0 * pow(3.0 * PI * PI, 2.0 / 3.0) * pow(rho, 5.0 / 3.0) / 10.0) + rhosq / (72 * rho) + lap / 6.0;
	double Dh = Fc*pow(rho,(5.0/ 3.0));
	kineng = kineng / 2.0 - (rhosq) / rho / 8.0;
	return 1 / (1 + (kineng / Dh)*(kineng / Dh));

}

double analysisBatch::RDG(double x, double y, double z)
{
	vectorReset(x, y, z);
	getDRho();
	double rhosq = dxrho*dxrho + dyrho*dyrho + dzrho*dzrho;


	if (rhosq == 0 || rho == 0)
	{
		return 999;
	}
	if (rho > rhoCutoff)
		return 999;

	return 0.161620459673995 * sqrt(rhosq) / (pow(rho, 4.0/3));
	//0.161620459673995D0 =  1/(2*(3*pi^2)^(1/3))
}

double* analysisBatch::AKinEng(double x, double y, double z,double * output)
{
	//derivInit(x, y, z);
	//wfnDerv();
	//wfnval();
	//wfnsdv();

	double rhosq = dxrho * dxrho + dyrho * dyrho + dzrho * dzrho;
	double lap = elecHess[0] + elecHess[4] + elecHess[8];
	output[0] = (3.0 * pow(3.0 * PI * PI, 2.0 / 3.0) * pow(rho, 5.0 / 3.0) / 10.0) + rhosq / (72 * rho) + lap / 6.0;
	output[1] = lap / 4.0 - 2 * output[0];
	output[2] = output[0] + output[1];
	return output;
}

double analysisBatch::vauleAtPoint(point Point, void * other)
{
	return RDG(Point.x * res + offsetx, Point.y * res + offsety, Point.z * res + offsetz);
}

void analysisBatch::cull(double x, double y, double z)
{
	vectorReset(x, y, z);
	std::vector<double> CX, CY, CZ;
	std::vector<int> CT;
	int C = 0;
	for(int i = 0; i < centers; i++)
	{
		if(distanceFromCenter[i] < 225)
		{
			CX.push_back(centerXvaule[i]);
			CY.push_back(centerYvaule[i]);
			CZ.push_back(centerZvaule[i]);
			CT.push_back(centertypes[i]);
			C++;
		}
	}
	centers = C;
	centerXvaule = new double[C];
	centerYvaule = new double[C];
	centerZvaule = new double[C];
	centertypes = new int[C];
	std::copy(CX.begin(), CX.end(), centerXvaule);
	std::copy(CY.begin(), CY.end(), centerYvaule);
	std::copy(CZ.begin(), CZ.end(), centerZvaule);
	std::copy(CT.begin(), CT.end(), centertypes);
	isCulled = 1;
}

/*
   double analysisBatch::elf(double x, double y, double z)
   {
   double rho = 0, Rhox = 0, rhoy = 0, rhoz = 0, rhosq;
   double kineng = 0;
   double Fc = 2.871234000; //magic number do not touch
   for (size_t i = 0; i < nmo; ++i)
   {
   kineng += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionDZ[i] * moWavefuntionDZ[i]);
   rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
   Rhox += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDX[i];
   rhoy += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDy[i];
   rhoz += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDZ[i];
   }
   Rhox *= 2;
   rhoy *= 2;
   rhoz *= 2;
   double Dh = Fc*pow(rho,(5.0/ 3.0));
   kineng = kineng / 2.0 - (Rhox*Rhox + rhoy*rhoy + rhoz*rhoz) / rho / 8.0;
   return 1 / (1 + (kineng / Dh)*(kineng / Dh));
   }*/
