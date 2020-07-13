#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include "array_tools.hpp"

#include "dcd_r.hpp"

using namespace std;

void DifferenceCalculation(const float* x, const float* y, const float* z, int numFrames, int currentFrame)
{
    vector<float> xvec (numFrames);
    vector<float> zvec (numFrames);
    vector<float> yvec (numFrames);
    vector<float> xvectemp (numFrames);
    vector<float> zvectemp (numFrames);
    vector<float> yvectemp (numFrames);

    double xdiff = 0;
    double zdiff = 0;
    double ydiff = 0;

	for(int u = 0; u < numFrames; u++)
	        {
	            //cout << "starting array copy of: " << x[u] << endl;
	            xvec.at(u) = x[u];
	            yvec.at(u) = y[u];
	            zvec.at(u) = x[u];
	            //cout << "Data that was copied in: " << xvec.at(u) << endl;
	        }
	if(currentFrame == 0)
	        {
	            //const float *tempx = x;
	            //tempx = copy(x, tempx);
	            //copy(tempx, x+100, tempx.begin())
	            //xvectemp = xvec;
	            for (int k = 0; k < numFrames; k++)
	            {
	                cout << "x coordinates of atom: " << k << "  :  " << xvec[k] << endl;
	            }
	        }
	        else
	        {

	            for(int k = 0; k < numFrames; k++)
	            {
	                //cout <<"Starting the difference calculation" << endl;
	                cout << "x coordinates of atom:" << k << " :      " << xvec[k] << endl;
	                cout << "X cocrdinates from last frame: " << xvectemp[k] << endl;

	                xdiff = xvec[k] - xvectemp[k];
	                zdiff = zvec[k] - zvectemp[k];
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
	                cout << "Difference in Y from last frame:    " << ydiff << endl;
	                cout << "Difference in Z from last frame:    " << zdiff << endl;
	                cout << " " << endl;
	                //xvectemp = xvec;
	                //copy(tempx, x+100, tempx.begin())
	            }
	        }

	        for(int b = 0; b < numFrames; b++)
	        {
	            cout << "starting temp array copy of: " << x[b] << endl;
	            xvectemp.at(b) = x[b];
	            zvectemp.at(b) = z[b];
	            yvectemp.at(b) = z[b];
	            cout << "Temp data that was copied in: " << xvectemp.at(b) << endl;
	        }
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
	cout <<"Number of first atom for analysis:" << endl;
	cin >> atom1;
	cout << "Number of second atom for analysis:" << endl;
	cin >> atom2;

	cout << "Variables are as follows: " << file << " ### " << atom1 << " ### " << atom2 << " ### " << endl;

	const char * filename = file.c_str();
	//filename = file;
	//DCD_R dcdf("newmd.dcd");
    DCD_R dcdf(filename);


    // read the header and print it
    dcdf.read_header();
    dcdf.printHeader();
    int numFrames = dcdf.getNFILE();
    const float *x,*y,*z;

    
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

        DifferenceCalculation(x,y,z,numFrames,i);
        //Start moving from here and call the function from here

        //frame counter
        cout << "For frame " << i << endl;

        //final print of header for additional information
        dcdf.printHeader();
        
        /* ... */
        
    }
    
    return EXIT_SUCCESS;
}










