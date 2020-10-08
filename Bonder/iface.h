#pragma once


#include "stdafx.h"
#include "point.h"
#include "readwfn.h"

//here be dragons

//weird mathy stuff 
//this pile of math is used for calculating the intracate math behind wave funtions
class analysisBatch
{
public:
	analysisBatch(wfnData input);
	~analysisBatch();
	//called by the flood fill ajusted from the grid to floating space
	double vauleAtPoint(point Point, void *other);
	//sets up the x,y,z offsets
	void setUpBatch(double x, double y, double z, double res);
	double elf();
	//reduced density gradent
	double RDG(double x, double y, double z);
	//arkinsoms reduced energy gradent
	double* AKinEng(double x, double y, double z,double * output);
	//reduced energy gradent and the signed vaule of rho
	double RDG_rho(double x, double y, double z, double *srho);
	void cull(double x, double y, double z);

	//x,y,z of atom i
	double atomx(int i);
	double atomy(int i);
	double atomz(int i);
private:
	//clears the vectors
	void vectorReset(double x, double y, double z);
	double getSignOfSecondEiganVaule();
	void getRho();
	void getDRho();
	void getDDRho();
	void getHesRho();

	double *__restrict__ dx, *__restrict__ dy, *__restrict__ dz;
	double *__restrict__ elecHess;
	double  rho,dxrho,dyrho,dzrho,ddxrho,ddyrho,ddzrho;


	int isCulled =0;
	double *__restrict__ centerXvaule, *__restrict__ centerYvaule, *__restrict__ centerZvaule,*__restrict__ distanceFromCenter;
	int *__restrict__ centertypes;
	//magic numbers
	double *__restrict__ czeta1,*__restrict__ czeta2,*__restrict__ czeta3,*__restrict__ czzeta1,*__restrict__ czzeta2,*__restrict__ czzeta3;
	double offsetx, offsety, offsetz, res;
	int centers;
};

static const int maxAtom = 18; // this code can only go up to argon
const double c1[maxAtom] = {0.2815, 2.437, 11.84, 31.34, 67.82, 120.2, 190.9,289.5,  406.3, 561.3, 760.8, 1016, 1319, 1658,2042, 2501, 3024, 3625};
const double c2[maxAtom] = {0,0, 0.06332, 0.3694, 0.8527, 1.172, 2.247,2.879, 3.049,6.984,22.42,37.17,57.95, 87.16,115.7, 158.0,   205.5,  260.0};
const double c3[maxAtom] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.06358, 0.3331, 0.8878, 0.7888, 1.465, 2.170,3.369, 5.211};
const double zeta1[maxAtom] = {0.5288, 0.3379, 0.1912, 0.1390, 0.1059, 0.0884,0.0767, 0.0669, 0.0608, 0.0549, 0.0496, 0.0449,0.0411, 0.0382, 0.0358, 0.0335, 0.0315, 0.0296};
const double zeta2[maxAtom] = {1, 1, 0.9992, 0.6945, 0.5300, 0.5480,0.4532, 0.3974, 0.3994, 0.3447, 0.2511, 0.2150,0.1874, 0.1654, 0.1509, 0.1369, 0.1259, 0.1168};
const double zeta3[maxAtom] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1.0236, 0.7753, 0.5962, 0.6995, 0.5851, 0.5149,0.4974, 0.4412};

