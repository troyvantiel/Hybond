#pragma once

#include "stdafx.h"
#include "LLI.h"
#include "iface.h"
#include "point.h"


struct gridPoint
{
	
	LLE* edges;
	LLi* internalPoints;	
	
	gridPoint()
	{
		edges = new LLE;
		internalPoints = new LLi;
	}

	void deletePoint()
	{
		delete edges;
		delete internalPoints;
	}
};

struct grid
{

	gridPoint *data;
	int x, y;
	~grid()
	{
		for (int i = 0; i < y; i++)
		{
			for (int j = 0; j < x; j++)
			{
				data[j * x +i].deletePoint();
			}
		}
		delete [] data;
	}
		
};

//runs the flood fil algorytem
grid fill(int x, int y, int z, void * other, int Xsize, int Ysize, double cutOff, bool *sucsess, analysisBatch* batch);
gridPoint* getPoint(grid* Grid, int x, int y);
