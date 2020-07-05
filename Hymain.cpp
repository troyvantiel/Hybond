#include <cstdlib>
#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

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
    
    const float *x,*y,*z;
    
    // in this loop the coordinates are read frame by frame
    for(int i=0;i< 500;i++)
    {
        cout<< "Getting Frame: " << i << endl;
        dcdf.read_oneFrame();
        cout<< "Finished Getting Frame: " << i << endl;
        
        /* your code goes here */

        
        //Getting x,y and z Co-ordinates and storing them in an array
        x=dcdf.getX();
        y=dcdf.getY();
        z=dcdf.getZ();

        if(i==0)
        {
            const float *tempx = null;
            for (int k = 0; k < 5; k++)
            {
                cout << "x coordinates of atom: " << k << "  :  " << x[k] << endl;
            }
        }
        else
        {
            tempx = x
            for(int k =0; k < 5; k++)
            {
                cout << "x coordinates of atom: " << k << "  :  " << x[k] << endl;
            }
        }
        cout << "For frame " << i << endl;

        //For loop that prints 5 x coords



        
        /* ... */
        
    }
    
    return EXIT_SUCCESS;
}