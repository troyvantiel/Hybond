#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <vector>

#include "array_tools.hpp"
#include "dcd_r.hpp"

using namespace std;


ofstream frameFile;
ofstream Diffoutput;


struct Coord
{
	double x;
	double y;
	double z;
};

void OutputDifferencetoFile(double xdiff, double ydiff, double zdiff, int frameCount)
{
	Diffoutput << "Difference in Frames:" << frameCount-1 << " and: " << frameCount << endl;
	Diffoutput << xdiff << ydiff << zdiff << endl;

}

void OutputFrametoFile(vector<Coord> atom, int frameCount, int numAtoms)
{

	 //all writing needs to happen after the open and before the close.

	frameFile << "Current Frame Printed: " << frameCount << endl;
	for(int i = 0; i < numAtoms; i++)
	{
		frameFile << atom[i].x << " " << atom[i].y << " " << atom[i].z << endl; //output coords of atoms in current frame to the file
	}

	//cout << "Frame outputted to file" << endl;
	//frameFile.close();
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
		for (int j = 1; j < numFrames; j++)
		{
			//cout << "Starting the difference calculation" << endl;
			//cout << "Current atoms Coordinates: X:" << atoms[j].x << "Y:" << atoms[j].y << "Z:" << atoms[j].z << endl;
			//cout << "Last atoms Coordinates: X:" << prevatom[j].x << "Y:" << prevatom[j].y << "Z:" << prevatom[j].z << endl;
			xdiff = atoms[j].x - prevatom[j].x;
			ydiff = atoms[j].y - prevatom[j].y;
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

	/*
    vector<float> xvec (numFrames);
    vector<float> zvec (numFrames);
    vector<float> yvec (numFrames);		//set up vectors to store all the data points from the frames
    vector<float> xvectemp (numFrames);
    vector<float> zvectemp (numFrames);
    vector<float> yvectemp (numFrames);

    double xdiff = 0;
    double zdiff = 0; //Set up the varibales to store the difference in the x,y,z coordinates
    double ydiff = 0;

    for(int u = 0; u < numFrames; u++)
    	{
	    	cout << "starting array copy of: " << x[u] << endl;
	        xvec.at(u) = x[u];
	        yvec.at(u) = y[u]; //array copying from const float* to vector<float> for ease of use
	        zvec.at(u) = x[u];
	        cout << "Data that was copied in: " << xvec.at(u) << endl;
	    }
	if(currentFrame == 0) //if the current frame is 0 no comparison can be made so alternate methods need to be taken
		{
			cout << "Current Frame is 0" << endl;
	    	for (int k = 0; k < numFrames; k++) //looping through all the data
	        	{
	                cout << "x coordinates of atom: " << k << "  :  " << xvec[k] << endl;
	            }
	    }
	else
		{
	    	for(int k = 0; k < numFrames; k++) //looping through the rest of the data
	        	{
	                //cout <<"Starting the difference calculation" << endl;
	                cout << "x coordinates of atom:" << k << " :      " << xvec[k] << endl;
	                cout << "X cocrdinates from last frame: " << xvectemp[k] << endl;
	                cout << "Y cocrdinates from last frame: " << yvectemp[k] << endl; //displaying of information
	                cout << "Z cocrdinates from last frame: " << zvectemp[k] << endl;


	                xdiff = xvec[k] - xvectemp[k];
	                zdiff = zvec[k] - zvectemp[k]; //calculation of the difference between the current frame and the last frame
	                ydiff = yvec[k] - yvectemp[k];

	                if(xdiff < 0)
	                {
	                	xdiff = xdiff * -1;//keep all differences positive values
	                }
	                if(zdiff < 0)
	                {
	                	zdiff = zdiff * -1;//keep all differences positive values
	                }
	                if(ydiff < 0)
	                {
	                	ydiff = ydiff * -1;//keep all differences positive values
	                }


	                cout << "Difference in X from last frame:    " << xdiff << endl;
	                cout << "Difference in Y from last frame:    " << ydiff << endl; //display all the differences
	                cout << "Difference in Z from last frame:    " << zdiff << endl;
	                cout << " " << endl;
	            }
	    	}

	        for(int b = 0; b < numFrames; b++)
	        {
	            cout << "starting temp array copy of: " << x[b] << endl;
	            xvectemp.at(b) = xvec[b];
	            zvectemp.at(b) = zvec[b]; //copy the current frame into the last frame arrays
	            yvectemp.at(b) = yvec[b];
	            cout << "Temp data that was copied in: " << xvectemp.at(b) << endl;
	        }
	        */
}

int main(int argc, char* argv[])
{                
    // instance of a new object DCD_R attached to a dcd file 
    string file = " ";
	int atom1 = 0;
	int atom2 = 0;

	//Get the user input for the dcd file and the atoms that want exploring
	cout << "Name of .dcd file you want to process: " <<endl;
	cin >> file;

	const char * filename = file.c_str(); //convert the filename into a string
	//filename = file;
	//DCD_R dcdf("newmd.dcd");
    DCD_R dcdf(filename); //read the dcd file with the filename that was given


    // read the header and print it
    dcdf.read_header();
    dcdf.printHeader();
    int numFrames = dcdf.getNFILE(); //get the number of frames from the header to read in
    int nAtom = dcdf.getNATOM();

	cout << "Number of atoms in the system: " << nAtom << endl;
	cout <<"Number of first atom for analysis:" << endl;
	cin >> atom1;
	cout << "Number of second atom for analysis:" << endl;
	cin >> atom2;

	cout << "Variables are as follows: " << file << " ### " << atom1 << " ### " << atom2 << " ### " << endl;


    const float *x,*y,*z; //make the const float varibles to store the coordinates.
    Coord atom;
    vector<Coord> atomsvec (numFrames);
    vector<Coord> lastvec (numFrames);

    frameFile.open ("Frames.txt"); // open file to be used by the frame output
	Diffoutput.open("DiffOut.txt");
    // in this loop the coordinates are read frame by frame
    for(int i=0; i < numFrames; i++)
    {


        //cout<< "Getting Frame: " << i << endl;
        dcdf.read_oneFrame();
        //cout<< "Finished Getting Frame: " << i << endl;
        
        /* your code goes here */

        
        //Getting x,y and z Co-ordinates and storing them in an array
        x = dcdf.getX();
        y = dcdf.getY();
        z = dcdf.getZ();

        //change the x,y,z coordinates into an atom struct that holds all that data


        for(int k = 0; k < nAtom; k++)
        {
        	atom.x = x[k];
        	atom.y = y[k];
        	atom.z = z[k];

        	atomsvec.at(k) = atom;
        }


        //Start moving from here and call the function from here
        lastvec = DifferenceCalculation(atomsvec, numFrames,i,lastvec);



        OutputFrametoFile(atomsvec, i, nAtom);






        //frame counter
        cout << "Finished frame: " << i << endl;

        //final print of header for additional information
        dcdf.printHeader();
        
        /* ... */
        
    }
    frameFile.close(); // close the frame output file to stop any bugs
	Diffoutput.close();
    return EXIT_SUCCESS;
}










