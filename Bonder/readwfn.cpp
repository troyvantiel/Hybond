#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "readwfn.h"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>
#include <string>

using namespace std;

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (!item.empty())
			elems.push_back(item);
	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

//decode fortran doubles
double DFD(string input)
{;
	std::vector<std::string> seglist = split(input,'D');
	double mines = stod(seglist[0]);
	double power = stoi(seglist[1]);
	return mines * pow(10, power);

}


int selectAtom(std::string name,int line)
{
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	if (name == "h") return 0;
	if (name == "he") return 1;
	if (name == "li") return 2;
	if (name == "be") return 3;
	if (name == "b") return 4;
	if (name == "c") return 5;
	if (name == "n") return 6;
	if (name == "o") return 7;
	if (name == "f") return 8;
	if (name == "ne") return 9;
	if (name == "na") return 10;
	if (name == "mg") return 11;
	if (name == "al") return 12;
	if (name == "si") return 13;
	if (name == "p") return 14;
	if (name == "s") return 15;
	if (name == "cl") return 16;
	if (name == "ar") return 17;
	cout << "atom not known, only atoms up to Ar are suported, error on line:" << line+2<< std::endl;
	cout << "atom read is: " << name << std::endl;
	exit(1);
}



wfnData* readFile(string file)
{
	
	wfnData* output = new wfnData;
	ifstream inputFile(file);
	
	if (!inputFile)
	{
		printf("no input\n");
		exit(1);
	}

	string line;
	std::getline(inputFile, line);
	output->nuc = stoi(line);
	std::getline(inputFile, line);
	output->x = new double[output->nuc];
	output->y = new double[output->nuc];
	output->z = new double[output->nuc];
	output->type = new int[output->nuc];
	output->name = new string[output->nuc];
	std::vector<std::string> tokens;
	for (int i = 0;i < output->nuc; i++)
	{
		std::getline(inputFile, line);
		tokens = split(line,' ');
		//conver to au
		output->x[i] = stod(tokens[1]) / 0.52917721092;
		output->y[i] = stod(tokens[2]) / 0.52917721092;
		output->z[i] = stod(tokens[3]) / 0.52917721092;
		output->type[i] = selectAtom(tokens[0],i);
		output->name[i] = tokens[0];
	}
	return output;
}	
