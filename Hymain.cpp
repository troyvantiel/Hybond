#include <cstdlib>
#include <iostream>

#include "array_tools.hpp"

#include "dcd_r.hpp"

using namespace std;

int main(int argc, char* argv[])
{                
    // instance of a new object DCD_R attached to a dcd file 
    DCD_R dcdf("prod.dcd");
    
    // read the header and print it
    dcdf.read_header();
    dcdf.printHeader();
    
    const float *x,*y,*z;
    
    // in this loop the coordinates are read frame by frame
    for(int i=0;i<dcdf.getNFILE();i++)
    {
        dcdf.read_oneFrame();
        
        /* your code goes here */
        
        x=dcdf.getX();
        y=dcdf.getY();
        z=dcdf.getZ();
        
        /* ... */
        
    }
    
    return EXIT_SUCCESS;
}