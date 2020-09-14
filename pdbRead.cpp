
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
	//cout << "Outputting the end of the vector" << endl;
	//cout << "size of infovec" << infovec.size() << endl;
	for(int i = 0; infovec.size(); i++)
	{
		if(infovec.back() == "Placeholder")
		{
			atomtypes.push_back("H");
			//cout << "Pushed a Hydrogen atom to the types due to placeholder found" << endl;
			break;
		}
		else
		{
			atomtypes.push_back(infovec.back());
			//cout << "Pushed non placeholder to the back of the atom types" << endl;
			break;
		}
	}
}




int main(int argc, char* argv[])
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
			int notemptycount =0;
			int whilecounter = 0;
			//cout << line << "\n";
			istringstream iss(line);
			//cout << line << endl;
			string s;
			while(getline(iss,s,' '))
			{
				//cout << "**WC: " << whilecounter << "**";
				if(whilecounter == 0)
				{
					if(s != "ATOM")
					{
						break;
					}
				}
				whilecounter++;
				count++;
				if(!s.empty())
				{
					notemptycount++;
					//cout << "**NEC:" << notemptycount << "**";
					if(notemptycount == 3)
					{
						cout << "Atom Types: " << s.front() << endl;

					}
				}
				//cout << count << s << endl;
				//cout << "Data at each split:" << s << endl;
				//cout << "Data at position 3 in the file" << s[2] << endl;
			}
			PasstoAtomData(infovec);
		}
		newfile.close();
	}
	for(int i =0; i < 20; i++ ) //needs to use atomtypes.size()
	{
		cout << i << atomtypes.at(i) << endl;
	}
	return 0;
}



