#pragma once
#include "readwfn.h"
#include "iface.h"

//this class takes an input point and anayses it
class  analysis
{
public:
	 analysis();
	~ analysis();
	void anilizePoint(int x, int y, int z, void * other, int Xsize, int Ysize, double cutOff, bool *sucsess, wfnData *data, std::string outputFile, analysisBatch* batch,int makeCube);
	void setUpAnalysisBatch(double x, double y, double z, double resalution, analysisBatch* batch);
private:
	double offsetx, offsety, offsetz, res;
};


