#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <fstream>            //simple c++ libs
#include <algorithm>
#include <vector>
#include <sstream>
#include <istream>
#include "array_tools.hpp"	  //hpp files for the arrays that the frame xyz information is stored in
#include "dcd_r.hpp"		  //dcd reader hpp file which allows reading of .dcd files
#include "Compat.hpp"         //hpp file to make Bonder methods available




using namespace std;

vector<string> atomtypes (0);
vector<double> atomdist(0);
ofstream frameFile;				//global variables for file output and the types of the atoms from the pdb file
ofstream Diffoutput;

struct Coord //struct for the Coord object. Holds all the information about the atom that might be needed for bonder to run.
{
	string type;	//this information is not stored in the .dcd file and is required from the .pdb file
	int atnum;      //number of atom within the molecule (APOA1 has 92224 atoms so this is from 0 to 92223)
	double x;
	double y;
	double z;
};

vector<string> readpdb()    //simple file reader for the .pdb file (it is in text format)
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
			int count =0;
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
					atomtypes.push_back(s);
					count++;
				}
			}
		}
		newfile.close();
	}
	return atomtypes;
}

void OutputDifferencetoFile(double xdiff, double ydiff, double zdiff, int frameCount)
{
	//outputting the differences to a file
	Diffoutput << "Difference in Frames:" << frameCount-1 << " and: " << frameCount << endl;
	Diffoutput << xdiff << ydiff << zdiff << endl;
}
void OutputRawFrametoFile(int currFrame, int numAtoms, vector<Coord> frames)
{
	ofstream rawfile;
	string filename = "./Frames/Frame:";
	filename += std::to_string(currFrame);
	filename += ".txt";
	rawfile.open(filename);
	rawfile << numAtoms-1 << endl;
	for(int i = 0; i<numAtoms; i++)
	{
		//string temptype ="";
		//temptype = frames[i].type;
		rawfile << frames[i].type << " " << frames[i].x << " " << frames[i].y << " " << frames[i].z << endl;
	}
	rawfile.close();
}
void OutputFrametoFile(vector<Coord> atom, int frameCount, int numAtoms)
{
	 //all writing needs to happen after the open and before the close.
	//frameFile << "Current Frame Printed: " << frameCount << endl;
	frameFile << numAtoms << endl;
	for(int i = 0; i < numAtoms; i++)
	{
		if(atom[i].x != 0 && atom[i].y != 0 && atom[i].z) //stop zeros beinng written to file (vector is filled up with 0,0,0)
		{
			frameFile << atom[i].type << " " <<atom[i].atnum << " " << atom[i].x << " " << atom[i].y << " " << atom[i].z << endl; //output coords of atoms in current frame to the file
		}

	}
}
vector<Coord> DifferenceCalculation(vector<Coord> atoms, int numFrames, int currentFrame, vector<Coord> prevatom)
{
	double xdiff = 0;
	double zdiff = 0;
	double ydiff = 0;
	if (currentFrame == 0)
	{
		//cout << "Current Frame is 0" << endl;
		for (int k = 0; k < numFrames; k++)
		{
			//cout << "Current atoms Coordinates: X:" << atoms[k].x << "Y:" << atoms[k].y << "Z:" << atoms[k].z << endl;
		}
	}
	else
	{
		for (int j = 0; j+1 < numFrames; j++)
		{
			//cout << "Starting the difference calculation" << endl;
			//cout << "Current atoms Coordinates: X:" << atoms[j].x << "Y:" << atoms[j].y << "Z:" << atoms[j].z << endl;
			//cout << "Last atoms Coordinates: X:" << prevatom[j].x << "Y:" << prevatom[j].y << "Z:" << prevatom[j].z << endl;
			xdiff = atoms[j].x - prevatom[j].x;
			ydiff = atoms[j].y - prevatom[j].y; //calculating all the differences of the atoms from the previous frame
			zdiff = atoms[j].z - prevatom[j].z;

			if(xdiff < 0)
			{
				xdiff = xdiff * -1;//keep all differences positive values
			}
			else if(zdiff < 0)
			{
				zdiff = zdiff * -1;//keep all differences positive values
			}
			else if(ydiff < 0)
			{
				ydiff = ydiff * -1;//keep all differences positive values
			}
			OutputDifferencetoFile(xdiff,ydiff,zdiff,currentFrame);
            //cout << "Difference in X from last frame:    " << xdiff << endl;
            //cout << "Difference in Y from last frame:    " << ydiff << endl; //display all the differences
            //cout << "Difference in Z from last frame:    " << zdiff << endl;
            //cout << " " << endl;
		}
	}
	return atoms;
}

int main(int argc, char* argv[])
{
	try{
		// instance of a new object DCD_R attached to a dcd file
		string file = " ";
		int atom1 = 0;
		int atom2 = 0;
		char version;
		char diff;
		char skip;
		cout << "Do you want to skip processing a .dcd and .pdb file?  (Y/n)" << endl;
		cin >> skip;

		if(skip != 'Y')
		{
				//Get the user input for the dcd file and the atoms that want exploring
				cout << "Name of .dcd file you want to process: " <<endl;
				cin >> file;
				//cout << "Software used to create file: ('N' for NAMD or 'C' for CHARMM)" << endl;
				//cin >> version;
				cout << "Enable Difference Output: (Warning large file size) Y/n"<< endl;
				cin >> diff;

				//convert the filename into a string and read dcd file from that name
				const char * filename = file.c_str();
				DCD_R dcdf(filename);

				// read the header and print it
				dcdf.read_header();
				dcdf.printHeader();

				//int numFrames = dcdf.getNPRIV();

				//get the number of frames from the header to read in
				int numFrames = dcdf.getNFILE();
				int nAtom = dcdf.getNATOM();


				cout << "Number of atoms in the system: " << nAtom << endl;
				cout <<"Number of first atom for analysis:" << endl;
				cin >> atom1;												//select the atoms to be looked at in the analysis
				cout << "Number of second atom for analysis:" << endl;
				cin >> atom2;

				//prints the variables that were chosen by the user
				cout << "Variables are as follows: " << file << " ### " << atom1 << " ### " << atom2 << " ### " << endl;

				 //make the const float varibles to store the coordinates.
				const float *x,*y,*z;
				Coord atom;

				//make the vectors to store all the information
				vector<Coord> atomsvec (nAtom);
				vector<Coord> lastvec (nAtom);
				vector<Coord> refinedVec (nAtom);
				vector<string> pdbVec (nAtom);

				//Calls the pead pdb file method to get the atom tpyes for analysis
				pdbVec = readpdb();

				//Open files to be used by the frame and difference output
				frameFile.open ("Frames.txt");
				Diffoutput.open("DiffOut.txt");
				// in this loop the coordinates are read frame by frame
				for(int i=0; i < numFrames; i++)
				{

					//Reads frame one by one and processes it
					//cout<< "Getting Frame: " << i << endl;
					dcdf.read_oneFrame();
					//cout<< "Finished Getting Frame: " << i << endl;

					//Your Code Goes Here


					//Getting x,y and z Co-ordinates and storing them in an array
					//cout << "Getting the x y and z coords" << endl;
					x = dcdf.getX();
					y = dcdf.getY();
					z = dcdf.getZ();
					//cout << "Finished getting the x y z coords for the frame" << endl;
					//cout << "for loop will run this many times:" << nAtom << endl;
					//cout << "length of each of the arrays holding the coords" << endl;

					//change the x,y,z coordinates into an atom struct that holds all the data
					for(int k = 0; k < nAtom; k++)
					{
						atom.type = pdbVec[k];
						atom.atnum = k;
						atom.x = x[k];
						atom.y = y[k];
						atom.z = z[k];
						//adds atom to the vector
						atomsvec.at(k) = atom;
					}



					//calculate the distance between the selected atoms and the atoms in the molecule using Euclidean Distance

					//set up coordinate variable for each atom
					Coord firstAtom = atomsvec.at(atom1);
					Coord secondAtom = atomsvec.at(atom2);
					//offset for adding to the file
					int dataAdd = 2;
					double sqx = (firstAtom.x - secondAtom.x) * (firstAtom.x - secondAtom.x);
					double sqy = (firstAtom.y - secondAtom.y) * (firstAtom.y - secondAtom.y);
					double sqz = (firstAtom.z - secondAtom.z) * (firstAtom.z - secondAtom.z);

					double dist = sqrt((sqx + sqz + sqy));
					atomdist.push_back(dist);

					for(int l = 0; l < nAtom; l++)
					{
						//start writing method for the code in the loop
						// variables that can be passed to the method: Each Atom , first Atom , second Atom , vector to add to , vector of atoms, dataAdd variable , loopcount

						Coord tempAtom = atomsvec.at(l);
						double fsqx = (firstAtom.x - tempAtom.x) * (firstAtom.x - tempAtom.x);
						double fsqy = (firstAtom.y - tempAtom.y) * (firstAtom.y - tempAtom.y); //finding the squared difference of each x y z for both the selected atoms
						double fsqz = (firstAtom.z - tempAtom.z) * (firstAtom.z - tempAtom.z);

						double ssqx = (secondAtom.x - tempAtom.x) * (secondAtom.x - tempAtom.x);
						double ssqy = (secondAtom.y - tempAtom.y) * (secondAtom.y - tempAtom.y);
						double ssqz = (secondAtom.z - tempAtom.z) * (secondAtom.z - tempAtom.z);

						double firstDist = sqrt((fsqx) + (fsqy) + (fsqz)); //calculating the distance the atom is from the selected ones
						double secondDist = sqrt((ssqx) +(ssqy) + (ssqz));

						if(firstDist <= 5 || secondDist <= 5) //filtering out all the needed atoms
						{
							refinedVec.at(0) = atomsvec.at(atom1); //adding the original picked atoms before the rest of the considered atoms
							refinedVec.at(1) = atomsvec.at(atom2);
							refinedVec.at(dataAdd) = atomsvec.at(l); //adding the new atom that is also to be considered in the future calculations
							dataAdd ++;
						}
					}
					//cout << "after the euclidean distance calculation" << endl;
					if(diff == 'Y')
					{
						//Start moving from here and call the function from here
						lastvec = DifferenceCalculation(atomsvec, numFrames,i,lastvec);
					}


					//cout << "if seen the error lies within the outputing frame to file function" << endl;
					//outputting the frame to the file with the new atoms that have been filtered out by distance
					OutputRawFrametoFile(i, nAtom, atomsvec);
					OutputFrametoFile(refinedVec, i, nAtom);

					//timing start
					//Bonder();
					//timing end;



					//frame counter
					cout << "Finished frame: " << i << endl;

					//final print of header for additional information
					//dcdf.printHeader();
					//cout << "After the print header at the end" << endl;

					/* ... */

				}
				frameFile.close(); // close the frame output file to stop any bugs
				Diffoutput.close();
		}
		for(int i = 0; i < 250; i++)
		{
			std::string newfile = "Frames/Frame:";
			newfile += std::to_string(i);
			newfile += ".txt";
			bond(argc, argv, newfile, i);
		}
	}
	catch(std::bad_alloc& ba)
	{
		std::cerr << "Bad alloc Caught: " << ba.what() << endl;
	}
    return EXIT_SUCCESS;
}










