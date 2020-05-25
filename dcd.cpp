
#include <cstdlib>
#include <iostream>

#include "dcd.hpp"

using namespace std;

DCD::DCD()
{
}

/*
 * When writing a block of binary data to a file, Fortran adds 2 unsigned integer before and after the block : each one contains the total number of bytes of the block.
 * 
 * For example for a fortran array of real(8) of size 10, 80 bytes of a data have to be written : so 3 things are written to the binary file : 
 *  1) An unsigned integer (4 bytes) , its value is 80
 *  2) The fortran array (80 bytes)
 *  3) A second unsigned integer (4 bytes), same value of 80.
 * This corresponds to the fortran statement "write(file_descriptor) fortran_array"
 * 
 * Now when writing several things at the same tine, i.e. "write(file_descriptor) fortran_array,an_integer"
 * things are written as following:
 * 1) unsigned int size
 * 2) fortran array
 * 3) an_integer
 * 4) unsigned int size
 * The value of the 2 unsigned ints is 84 : 80 for the array, 4 for the integer(4).
 * 
 * The following method DCD::checkFortranIOerror check that the 2 unsigned ints have the same value: if not there was a probem when reading the binary fortran file.
 */
void DCD::checkFortranIOerror(const char file[], const int line, const unsigned int fortcheck1, const unsigned int fortcheck2) const
{
    if( fortcheck1 != fortcheck2 )
    {
        cout << "Error when reading data from dcd : quantities do not match." << endl;
        cout << "fortcheck1 = " << fortcheck1 <<  " and fortcheck2 = " << fortcheck2 << endl;
        cout << "in File " << file << " at Line " << line << endl;
        exit(EXIT_FAILURE);
    }
}


int DCD::getNFILE() const {
    return NFILE;
}

const float* DCD::getZ() const {
    return Z;
}

const float* DCD::getY() const {
    return Y;
}

const float* DCD::getX() const {
    return X;
}

const int* DCD::getFREEAT() const {
    return FREEAT;
}

int DCD::getLNFREAT() const {
    return LNFREAT;
}

int DCD::getNATOM() const {
    return NATOM;
}

int DCD::getCHARMV() const {
    return CHARMV;
}

int DCD::getQCRYS() const {
    return QCRYS;
}

int DCD::getDELTA4() const {
    return DELTA4;
}

int DCD::getFROZAT() const {
    return FROZAT;
}

int DCD::getNDEGF() const {
    return NDEGF;
}

int DCD::getNSTEP() const {
    return NSTEP;
}

int DCD::getNSAVC() const {
    return NSAVC;
}

int DCD::getNPRIV() const {
    return NPRIV;
}

const double* DCD::getPbc() const {
    return pbc;
}

const char* DCD::getTITLE() const {
    return TITLE;
}

int DCD::getNTITLE() const {
    return NTITLE;
}

const int* DCD::getICNTRL() const {
    return ICNTRL;
}

const char* DCD::getHDR() const {
    return HDR;
}

DCD::~DCD()
{
}