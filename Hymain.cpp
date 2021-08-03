#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <time.h>
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


//global variables
vector<string> atomtypes (0); //vector for the atom types from the pdb file
vector<double> atomdist(0);   //vector to store the distance between the target atoms for each frame
ofstream frameFile;			  //output stream for the frames
ofstream Diffoutput;		  //output stream for the difference output to detect all atoms in a certain radius of the target atoms
clock_t t1, t2;				  //clock variables to measure the time it takes to call each iteration of Bonder

struct Coord //struct for the Coord object. Holds all the information about the atom that might be needed for bonder to run.
{
	string type;	//this information is not stored in the .dcd file and is required from the .pdb file
	int atnum;      //number of atom within the molecule (APOA1 has 92224 atoms so this is from 0 to 92223)
	double x;
	double y;		//doubles for storing the x,y and z coordinates
	double z;
};

void OutputTime(int frames, vector<double> timing) //method to output the time taken to a file
{
	ofstream timeOut;							   //output stream for the timing file
	string filename = "Runtimes.csv";			   //name of the timing file
	timeOut.open(filename);                        //opening the time file
	for(int i =0; i < frames; i++)				   //for loop to output each time in the array
	{
		timeOut << i << "," <<timing.at(i) << endl; //adding the data to the file
	}


}

void bondAngle() //method signature to be called to calculate the angle of the hydrogen bond as this has been shown to effect the energy of the bond
{

}

void OutputDistancetoFile(vector<double> atomdist, int numFrames) //method to output the distance between the target atoms
{
	cout << "Outputting distances to a file" << endl;
	ofstream distOut; 											//output stream for the distance file

	bool atomDistFlag = true;
	int longDistCount = 0;
	string distInput = "";

	string filename = "";										//filename variable
	//cout << "Output filename for the atom distances" << endl;
	cin >> filename;											//getting filename from the user
	distOut.open(filename);										//opening the distance file
	for(int i =0; i < numFrames; i++)							//for loop to output each distance
	{
		if(atomDistFlag == true)
		{
			if(atomdist.at(i) > 5)
			{
				longDistCount++;
				if(longDistCount > 30)
				{
					cout << "30 long distance hydrogen bonds have been detected" << endl;
					cout << "Do you want to Continue? (y/n)" << endl;//ask the user if they want to continue knowing that the bond is not well defined
					cin >> distInput;
					if(distInput.compare("n"))
					{
						exit(0);
						//if the answer is no kill the program if it is yes carry on.
					}
					else
					{
						atomDistFlag = false;
					}

				}
			}
		}

		distOut << atomdist.at(i) << endl;						//adding data to the file
	}

}

vector<string> readpdb()    //simple file reader for the .pdb file (it is in text format)
{
	string filename;									//string for storing the name of the file
	cout << "Type the filename of pdb file" << endl;
	cin >> filename;									//getting user input of the file
	fstream newfile;									//creating the stream to read the file
	newfile.open(filename, ios::in);					//opening the file
	if(newfile.is_open())								//if loop to check if the file is opened
	{
		string line;									//string for storing each line of data
		while(getline(newfile, line))					//while loop to get the line from the file and store it in the variable.
		{
			int count =0;								//counter for the while loop
			istringstream iss(line);					//converts the line of the file to a string stream
			string s;									//string variable to store each value of the string stream
			while(getline(iss,s,' '))					//while loop to read each token of the string stream
			{
				while(count == 0)						//while loop to check if the count is not 0
				{
					if(s.length() > 1)					//checks for the column that has atom type in it
					{
						break;							//breaks the while loop early when the string is too long to be an atom type
					}
					atomtypes.push_back(s);				//push the string containing the atom type to the back of the atomtypes vector
					count++;							//count variable to stop the while loop for that line of data
				}
			}
		}
		newfile.close();								//closing the file to stop I/O errors
	}
	return atomtypes;
}

void OutputDifferencetoFile(double xdiff, double ydiff, double zdiff, int frameCount)
{
	//outputting the differences to a file
	Diffoutput << "Difference in Frames:" << frameCount-1 << " and: " << frameCount << endl; //outputting the difference to a file
	Diffoutput << xdiff << ydiff << zdiff << endl;											 //outputting the actual difference values
}
void OutputRawFrametoFile(int currFrame, int numAtoms, vector<Coord> frames) //method for outputting raw frames to files
{
	ofstream rawfile;						//output stream
	string filename = "./Frames/Frame:";	//filename for file
	filename += std::to_string(currFrame);	//adding the current frame number to file name
	filename += ".txt";						//adding the last part of the filename
	rawfile.open(filename);					//open the raw frame file
	rawfile << numAtoms-1 << endl;			//output the number of atoms with one index taken off to stop out of bounds
	for(int i = 0; i<numAtoms; i++)			//for loop to go through each frame
	{
		//string temptype ="";
		//temptype = frames[i].type;
		rawfile << frames[i].type << " " << frames[i].x << " " << frames[i].y << " " << frames[i].z << endl; //output each data point with the correct formatting
	}
	rawfile.close(); //close the file to prevent I/O output
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
vector<Coord> DifferenceCalculation(vector<Coord> atoms, int numFrames, int currentFrame, vector<Coord> prevatom) //output for the difference of atoms arounf the target ones
{
	double xdiff = 0;
	double zdiff = 0; //doubles to store the difference in the x,y and z values
	double ydiff = 0;
	if (currentFrame == 0) //if the current frame is the first one skip beacuse there is no difference calculation for the first frame
	{
		return atoms; //if the current frame is 0 there is only one frame that has been read in by the program and so nothing needs to be calculated.
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
		string file = " ";
		int atom1 = 0;
		int atom2 = 0;
		char version;           //variables to store various data and flags for the program
		char diff;
		char skip;
		double time =0;
		cout << "Do you want to skip processing a .dcd and .pdb file?  (Y/n)" << endl; //asking the user if they want to skip the processing of a .dcd as this only needs to happen once for each simulation
		cin >> skip; //getting the user input

		if(skip != 'Y') //checking the user input to determine if .dcd processing is needed
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
					dcdf.read_oneFrame();

					//######### Any code that will manipulate each frame that comes from the .dcd file can be written here ##########


					//Getting x,y and z Co-ordinates and storing them in an array
					x = dcdf.getX();
					y = dcdf.getY();
					z = dcdf.getZ();

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
					//cout << firstAtom.x << "  " << firstAtom.y << "  " << firstAtom.z << endl;
					//cout << secondAtom.x << "  " <<secondAtom.y << "  " << secondAtom.z << endl;

					//set up for calculating the Euclidean distance between the two target atoms
					double sqx = (firstAtom.x - secondAtom.x) * (firstAtom.x - secondAtom.x);
					double sqy = (firstAtom.y - secondAtom.y) * (firstAtom.y - secondAtom.y);
					double sqz = (firstAtom.z - secondAtom.z) * (firstAtom.z - secondAtom.z);
					//calculating the Euclidean distance for the target atoms for each frame
					double dist = sqrt((sqx + sqz + sqy));
					//push the distance to the back of an array
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
					if(diff == 'Y') //checking the value of the "diff" flag
					{
						//Start moving from here and call the function from here
						lastvec = DifferenceCalculation(atomsvec, numFrames,i,lastvec);
					}
					//outputting the frame to the file with the new atoms that have been filtered out by distance
					OutputRawFrametoFile(i, nAtom, atomsvec);
					OutputFrametoFile(refinedVec, i, nAtom);
					//frame counter and output for user feedback
					if (i % 10 == 0)
					{
						cout << "Finished frame: " << i << endl;
					}
				}
				//call to output the distance of the two target atoms
				OutputDistancetoFile(atomdist, numFrames);
				frameFile.close(); // close the frame output file to stop any bugs
				Diffoutput.close();
		}
		//vector to store the runtimes of each Bonder call
		vector<double> timing(0);
		//for loop to run foreach frame of the file
		for(int i = 0; i < 250; i++)
		{
			cout << "Current Frame Being Processed: " << i << endl;
			std::string newfile = "Frames/Frame:"; //string for the filename
			newfile += std::to_string(i);		   //adding the number and the back end of the filename
			newfile += ".txt";
			t1 = clock();						//take the current program time just before the call to bonder
			bond(argc, argv, newfile, i);
			t2 = clock();						//take the current program time just after the call to bonder
			time = (t2 - t1)/(CLOCKS_PER_SEC/(double) 1000.0);	//calculate the time in milliseconds
			timing.push_back(time);								//push the data onto a vector for outputting to a file
			cout << "Runtime of the Bonder call for the Frame: " << time << endl;		//print statement to show how long bonder takes each time
		}
		OutputTime(250, timing); //outputting the timings to a file
	}
	catch(std::bad_alloc& ba) //catch statement for any bad I/O errors
	{
		std::cerr << "Bad alloc Caught: " << ba.what() << endl; //print the stack trace
	}
    return EXIT_SUCCESS; //return statement
}










