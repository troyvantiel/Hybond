#include "dcd.hpp"

#ifndef DCD_R_HPP
#define	DCD_R_HPP

class DCD_R : public DCD
{

private:
    //no private attributes
    //private methods
    void alloc();
    
public:
    
    // no public attributes
    // public methods
    DCD_R(const char filename[]); //constructor
    //DCD_R(string filename);
    
    void read_header(char check);
    void read_oneFrame(char check);
    void printHeader() const;
        
    ~DCD_R();

};

#endif	/* DCD_R_HPP */
