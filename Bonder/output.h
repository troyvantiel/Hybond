#pragma once
#include "readwfn.h"
#include "iface.h"
void outputCube(double minx, double miny, double minz, double maxx, double maxy, double maxz, double res, std::string file, wfnData inputData,double cutoff,analysisBatch* batch,int makeCube);
