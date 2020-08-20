#include "pdbReader.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <list>

using namespace std;
vector<string> atomtypes (0);

void PasstoAtomData(vector<string> infovec)
{
	cout << "Outputting the end of the vector" << endl;
	cout << "size of infovec" << infovec.size() << endl;
	/*for(int i  =0; infovec.size(); i++)
	{
		if(infovec.back() == "Placeholder")
		{
			infovec.push_back("H");
			cout << "Pushed a Hydrogen atom to the types due to placeholder found" << endl;
		}
		else
		{
			atomtypes.push_back(infovec.back());
			cout << "Pushed non placeholder to the back of the atom types" << endl;
		}
	}*/
}




int main()
{
	string filename;
	vector<string> infovec (1);
	cout << "Type the filename of pdb file" << endl;
	cin >> filename;
	fstream newfile;
	newfile.open(filename, ios::in);
	if(newfile.is_open())
	{
		string line;
		while(getline(newfile, line))
		{
			int count =0;
			//cout << line << "\n";
			istringstream iss(line);
			string s;
			while(getline(iss,s,' '))
			{
				if(s.empty() == true)
				{
					s = "Placeholder";
				}
				//cout << count << s << endl;
				infovec.push_back(s);
				count++;
			}
			PasstoAtomData(infovec);
		}
		newfile.close();
	}
	for(int i =0; i < infovec.size(); i++ )
	{

	}
}


