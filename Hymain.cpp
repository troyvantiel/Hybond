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

int main(int argc, char* argv[])
{                
    // instance of a new object DCD_R attached to a dcd file 
    DCD_R dcdf("newmd.dcd");
    
    // read the header and print it
    dcdf.read_header();
    dcdf.printHeader();
    
    const float *x,*y,*z, *tempx;
    vector<float> xvec (200);
    vector<float> xvectemp (200);
    
    // in this loop the coordinates are read frame by frame
    for(int i=0;i < 4;i++)
    {
        //cout<< "Getting Frame: " << i << endl;
        dcdf.read_oneFrame();
        //cout<< "Finished Getting Frame: " << i << endl;
        
        /* your code goes here */

        
        //Getting x,y and z Co-ordinates and storing them in an array
        x = dcdf.getX();


        y=dcdf.getY();
        z=dcdf.getZ();


        // for loop to copy data into a vector for ease of management
        for(int u =0; u<200; u++)
        {
            //cout << "starting array copy of: " << x[u] << endl;
            xvec.at(u) = x[u];
            //cout << "Data that was copied in: " << xvec.at(u) << endl;
        }


        //if statement to catch the first frame (no comparison)
        if(i==0)
        {
            //const float *tempx = x;
            //tempx = copy(x, tempx);
            //copy(tempx, x+100, tempx.begin())
            //xvectemp = xvec;
            for (int k = 0; k < 5; k++)
            {
                cout << "x coordinates of atom: " << k << "  :  " << xvec[k] << endl;
            }
        }
        else
        {

            for(int k = 0; k < 5; k++)
            {
                //cout <<"Starting the difference calculation" << endl;
                cout << "x coordinates of atom:" << k << " :      " << xvec[k] << endl;
                cout << "X cocrdinates from last frame: " << xvectemp[k] << endl;
                double diff = xvec[k] - xvectemp[k];
                if(diff<0)
                {
                	diff = diff *-1;//keep all differences positive values
                }
                cout << "Difference from last frame:    " << diff << endl;
                cout << " " << endl;
                //xvectemp = xvec;
                //copy(tempx, x+100, tempx.begin())
            }
        }

        for(int b = 0; b<200; b++)
        {
            //cout << "starting temp array copy of: " << x[b] << endl;
            xvectemp.at(b) = x[b];
            //cout << "Temp data that was copied in: " << xvectemp.at(b) << endl;
        }
        cout << "For frame " << i << endl;

        //For loop that prints 5 x coords


        dcdf.printHeader();
        
        /* ... */
        
    }
    
    return EXIT_SUCCESS;
}
