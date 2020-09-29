
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <list>

using namespace std;
vector<string> atomtypes (0);
//Note: This approach doesnt give enough information for determining
//Note:	the following atoms: He,Ne,Na,Be,Si,Cl
vector<string> elementcheck{"H","C","O","N","F","B","P","S","NA","MG","AL","SI","CL","LI","BE","NE","AR","HE"};

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
				//cout << line << "\n";
				istringstream iss(line);
				string s;
				while(getline(iss,s,' '))
				{
					while(count == 0)
					{
						if(s.length() > 1)
						{
							break;
						}
						cout << count << s << endl;
						atomtypes.push_back(s);
						count++;
					}
				}

			}
			newfile.close();
		}
		for(int i =0; i < atomtypes.size(); i++ ) //needs to use atomtypes.size()
		{
			cout << "Vector *atomtypes* at position: " << i << " -- "<< atomtypes.at(i) << endl;
		}
		return 0;
}



