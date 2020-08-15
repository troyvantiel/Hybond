#include "pdbReader.hpp"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
	string filename;
	cout << "Type the filename of pdb file" << endl;
	cin >> filename;
	fstream newfile;
	newfile.open(filename, ios::in);
	if(newfile.is_open())
	{
		string line;
		while(getline(newfile, line))
		{
			cout << line << "\n";
		}
		newfile.close();
	}
}


