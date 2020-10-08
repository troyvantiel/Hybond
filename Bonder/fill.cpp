#include "stdafx.h"
#include "fill.h"
#include <queue>


void addToGrid(gridPoint* EditingPoint,point Point,void* other,double cutOff, analysisBatch* batch);
void getAllNeibors(point currentPoint,point *neibors);
bool checkIfUsed(gridPoint* editingPoint, point Point);
bool checkIfEdge(point Point, void *other, double cutoff, analysisBatch* batch);

void fillOnce(std::queue<point> *toProcess,grid *data,double cutOff, void* other,analysisBatch* batch,int Xsize,int Ysize)
{
	point currentPoint;
	currentPoint = (*toProcess).front();
	(*toProcess).pop();

	gridPoint* EditingPoint = getPoint(data, currentPoint.x, currentPoint.y);

	if (currentPoint.x == ((*data).x / 2) -3 || currentPoint.x == -((*data).x / 2) +3 || currentPoint.y == (*data).y / 2 -3 || currentPoint.y == -((*data).y / 2) +3)
		return;
	if (checkIfUsed(EditingPoint, currentPoint))
		return;
	if ((*batch).vauleAtPoint(currentPoint, 0) >= cutOff)
		return;
	if (!checkIfEdge(currentPoint, other, cutOff,batch))
		return;
	addToGrid(EditingPoint, currentPoint, other, cutOff,batch);


	//printf("%d %d %d \n",currentPoint.x,currentPoint.y,currentPoint.z);
	point newPoints[26];
	getAllNeibors(currentPoint,newPoints);
	for (size_t i = 0; i < 26; i++)
	{
		if (checkIfUsed(getPoint(data, newPoints[i].x, newPoints[i].y), currentPoint))
			continue;

		(*toProcess).push(newPoints[i]);
	}
}


grid fill(int x, int y, int z, void * other, int Xsize, int Ysize, double cutOff,bool *sucsess,analysisBatch* batch)
{
	grid data;
	data.x = Xsize;
	data.y = Ysize;
	data.data = new gridPoint[Xsize*Ysize]();
	point currentPoint;
	currentPoint.x = x;
	currentPoint.y = y;
	currentPoint.z = z;

	if ((*batch).vauleAtPoint(currentPoint, other) >= cutOff)
	{
		*sucsess = false;
		return data;
	}
	
	while ((*batch).vauleAtPoint(currentPoint, other) < cutOff)
	{
		currentPoint.z++;
	}
	currentPoint.z--;
	std::queue<point> toProcess;
	toProcess.push(currentPoint);
	int l = 0;
	while (!(toProcess.empty()))
	{
		
		fillOnce(&toProcess, &data, cutOff, other,batch,Xsize,Ysize);
	}

	*sucsess = true;
	return data;
}

gridPoint* getPoint(grid* Grid, int x, int y)
{
	int loc = (x + (*Grid).x / 2) + (*Grid).x * (y + (*Grid).y / 2);
	return (*Grid).data + loc;
}

//front is defined as z+1 back is z -1
void addToGrid(gridPoint* EditingPoint, point Point,void *other,double cutoff,analysisBatch* batch)
{
	Point.z++;
	bool front = (*batch).vauleAtPoint(Point, other) < cutoff;
	Point.z -= 2;
	bool back = (*batch).vauleAtPoint(Point,other) < cutoff;
	Point.z++;
	if (front && back)
	{
		(*(*EditingPoint).internalPoints).add(Point.z);

	}
	else
	{
		edgepoint toAdd;
		toAdd.Z = Point.z;
		if (front)
			toAdd.LR = 1;
		else if (back)
			toAdd.LR = 2;
		else
			toAdd.LR = 0;

		(*(*EditingPoint).edges).add(toAdd);
	}



}

void getAllNeibors(point currentPoint,point *neibors)
{
	point aPoint;
	int location;
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				location = i + 3 * j + 9 * k;
				if (location == 13)
					continue;
				if (location > 13)
				{
					location--;
				}
				aPoint.x = currentPoint.x + i - 1;
				aPoint.y = currentPoint.y + j - 1;
				aPoint.z = currentPoint.z + k - 1;
				neibors[location] = aPoint;
			}
		}
	}
}

point* getOrthoNeibors(point currentPoint,point *neibors)
{
	for (size_t i = 0; i < 6; i++)
	{
		neibors[i] = currentPoint;
	}
	neibors[0].x--;
	neibors[1].x++;
	neibors[2].y--;
	neibors[3].y++;
	neibors[4].z--;
	neibors[5].z++;
	return neibors;
}

bool checkIfUsed(gridPoint* editingPoint, point Point)
{
	gridPoint EditingPoint = *editingPoint;
		if ((*(EditingPoint.edges)).search(Point.z))
		{
			return true;
		}


		if ((*(EditingPoint.internalPoints)).search(Point.z))
		{
			return true;
		}

	return false;
}

bool checkIfEdge(point Point, void *other, double cutoff,analysisBatch* batch)
{
	point neibors[6];
	getOrthoNeibors(Point,neibors);
	for (size_t i = 0; i < 6; i++)
	{
		if ((*batch).vauleAtPoint(neibors[i], other) >= cutoff)
		{
			return true;
		}
			
	}
	return false;
}
