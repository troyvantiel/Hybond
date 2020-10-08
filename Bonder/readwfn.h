#pragma once
#include <string>
#include <iostream>

//this class reads in .wfn files
struct wfnData
{
	int nuc;
	
	//atoms
	double *x, *y, *z;
	int *type;
	std::string *name;

	//assinments

	//MO

};

wfnData* readFile(std::string file);
