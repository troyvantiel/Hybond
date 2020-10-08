#include "stdafx.h"
#include "fill.h"
#include <algorithm>
#include <vector>
#include "analize.h"
#include "output.h"
#include <mutex>
#include <list>
#include <iostream>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

boost::mutex mutexCube,mutexCenters;

struct center
{
	double x, y, z;
	center(double X, double Y, double Z)
	{
		x = X;
		y = Y;
		z = Z;
	}

	double disi(center cent)
	{
		return (x - cent.x) * (x - cent.x) + (y - cent.y) * (y - cent.y) + (z - cent.z) * (z - cent.z);
	}
	inline bool operator==(const center lhs);
};

inline bool center::operator==(const center lhs) { return (lhs.x == x) && (lhs.y == y) && (lhs.z == z); }

struct cube
{
	double maxX, maxY, maxZ, minX, minY, minZ;
	cube(double NX,double NY,double NZ, double XX, double XY, double XZ)
	{
		maxX = XX;
		maxY = XY;
		maxZ = XZ;
		minX = NX;
		minY = NY;
		minZ = NZ;
	}
};

std::vector<cube> *spacesUsed = new std::vector<cube>;
std::list<center> *centers = new std::list<center>;

//to sort edge points
bool edgeComp(edgepoint i, edgepoint j)
{
	return i.Z < j.Z;
}

void analysis::setUpAnalysisBatch(double x, double y, double z, double resalution, analysisBatch* batch)
{
	//x = res * ((int)(x / res));
	//y = res * ((int)(y / res));
	//z = res * ((int)(z / res));
	(*batch).setUpBatch(x, y, z, resalution);
	offsetx = x;
	offsety = y;
	offsetz = z;
	res = resalution;
}

void analysis::anilizePoint(int x, int y, int z, void * other, int Xsize, int Ysize, double cutOff, bool *sucsess,wfnData *data,std::string outputFile, analysisBatch* batch,int makeCube)
{
	//std::cout << makeCube << std::endl;
	bool looping = true;
	center thisPoint(offsetx, offsety, offsetz);
	while (looping)
	{
		mutexCube.lock();
		for (int i = 0; i < (*spacesUsed).size(); i++)
		{
			cube var = (*spacesUsed)[i];
			if (offsetx < var.maxX && offsetx > var.minX && offsety< var.maxY && offsety > var.minY && offsetz < var.maxZ && offsetz > var.minZ)
			{
				mutexCube.unlock();
				*sucsess = false;
				return;
			}
		}
		mutexCube.unlock();
		mutexCenters.lock();
		bool found = false;
		for (std::list<center>::iterator it = ((*centers).begin()); it != ((*centers).end()); ++it)
		{
			if ((*it).disi(thisPoint) < 5)
			{
				found = true;
				break;
			}
		}
		if (!found)
		{
			(*centers).push_back(thisPoint);
			looping = false;
			mutexCenters.unlock();
		}
		else
		{
			mutexCenters.unlock();
			boost::this_thread::sleep_for(boost::chrono::seconds(10));
		}
		
	}

	
	(*batch).cull(offsetx,offsety,offsetz);	
	printf("starting filling at %f %f %f\n", offsetx, offsety, offsetz);
	grid results = fill(x,y,z,other,Xsize,Ysize,cutOff,sucsess,batch);
	if (!*sucsess)
	{
		printf("fail\n");
		return;
	}
	printf("finished filling\n");
	int maxX=x, maxY=y, maxZ=0, minX=x, minY=y, minZ=0;
	int vol = 0;
	gridPoint currentPoint;
	for (int i = -Xsize/2; i < Xsize/2; i++)
	{
		for (int j = -Ysize/2; j < Ysize/2; j++)
		{
			currentPoint = *getPoint(&results, i, j);
			int size;
			edgepoint* edges = (*(currentPoint.edges)).dump(&size);
			if (size == 0)
				continue;
			if (i > maxX)
				maxX = i;
			if (i < minX)
				minX = i;

			if (j > maxY)
				maxY = j;
			if (j < minY)
				minY = j;
			
			std::sort(edges, edges + size, edgeComp);
			int lastZ = INT32_MAX;
			for (size_t k = 0; k < size; k++)
			{
				//printf("%d %d %d\n", i, j, edges[k].Z);
				if (edges[k].Z < minZ)
					minZ = edges[k].Z;
				if (edges[k].Z > maxZ)
					maxZ = edges[k].Z;

				switch (edges[k].LR)
				{
				case 0:
					vol++;
					break;
				case 1:
					lastZ = edges[k].Z;
					break;
				case 2:
					vol += edges[k].Z - lastZ + 1;
					lastZ = INT32_MAX;
				}
			}
			delete [] edges;

		}
	}

	minX -= 10;
	minY -= 10;
	minZ -= 10;

	maxX += 10;
	maxY += 10;
	maxZ += 10;	

	mutexCube.lock();
	for (int i = 0; i < (*spacesUsed).size(); i++)
	{
		cube var = (*spacesUsed)[i];
		if (offsetx < var.maxX && offsetx > var.minX && offsety< var.maxY && offsety > var.minY && offsetz < var.maxZ && offsetz > var.minZ)
		{
			mutexCube.unlock();
			*sucsess = false;
			return;
		}
	}
	(*spacesUsed).push_back(cube(minX* res + offsetx, minY* res + offsety, minZ* res + offsetz, maxX* res + offsetx, maxY* res + offsety, maxZ* res + offsetz));
	mutexCube.unlock();

	mutexCenters.lock();
	(*centers).remove(thisPoint);
	mutexCenters.unlock();

	printf("drawing grid at %f %f %f, %f %f %f\n", minX* res + offsetx, minY* res + offsety, minZ* res + offsetz, maxX* res + offsetx, maxY* res + offsety, maxZ* res + offsetz);
	printf("volume is  %f\n", vol*res*res*res);

	if (!(outputFile.empty()))
		outputCube(minX* res + offsetx, minY* res + offsety, minZ* res + offsetz, maxX* res + offsetx, maxY* res + offsety, maxZ* res + offsetz, res, outputFile, *data,cutOff,batch,makeCube);
	/*for (size_t i = 0; i < Xsize; i++)
	{
		for (size_t j = 0; j < Ysize; j++)
		{
			currentPoint = *getPoint(&results, i, j);
			delete (currentPoint.edges);
			delete (currentPoint.internalPoints);
		}
	}
	delete results.data;
	*/

}

analysis::analysis()
{
}

analysis::~analysis()
{
}

