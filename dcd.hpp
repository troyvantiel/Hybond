#ifndef DCD_HPP_INCLUDED
#define DCD_HPP_INCLUDED

#include <fstream>

class DCD
{
protected:
    //protected attributes
    std::fstream dcdf;  // for opening dcd file
    
    bool dcd_first_read; // at first read of the coordinates if there are some frozen atoms the number of coordinates to read is different than for other frames
    
    char HDR[4+1]; // is CORD if only coordinates or VEL if velocities included (not supported yet)
    int  ICNTRL[20];
    int  NTITLE; //how many "title lines" in dcd
    char *TITLE; //each "title line" is 80 char long.
    
    double pbc[6];  //a matrix of 6 real defining the periodic boundary conditions : only useful if QCRYS is not 0
    
    /*content of ICNTRL : non detailed ones are 0 */
    
    int NFILE;  //ICNTRL(1)  number of frames in this dcd
    int NPRIV;  //ICNTRL(2)  if restart, total number of frames before first print
    int NSAVC;  //ICNTRL(3)  frequency of writting dcd
    int NSTEP;  //ICNTRL(4)  number of steps ; note that NSTEP/NSAVC = NFILE
    int NDEGF;  //ICNTRL(8)  number of degrees of freedom
    int FROZAT; //ICNTRL(9) is NATOM - NumberOfFreeAtoms : it is the number of frozen (i.e. not moving atoms)
    int DELTA4; //ICNTRL(10) timestep in AKMA units but stored as a 32 bits integer !!!
    int QCRYS;  //ICNTRL(11) is 1 if CRYSTAL used.
    int CHARMV; //ICNTRL(20) is charmm version
    
    int NATOM; // Number of atoms
    
    int LNFREAT; // Number of free (moving) atoms.
    int *FREEAT; // Array storing indexes of moving atoms.
    
    //coordinates stored in simple precision (IMPORTANT)
    float *X;
    float *Y;
    float *Z;
    
    //protected methods
    virtual void alloc()=0;
    void checkFortranIOerror(const char file[], const int line, 
                             const unsigned int fortcheck1, const unsigned int fortcheck2) const;
public:
    DCD();
    
    int getNFILE() const;
    const float* getZ() const;
    const float* getY() const;
    const float* getX() const;
    const int* getFREEAT() const;
    int getLNFREAT() const;
    int getNATOM() const;
    int getCHARMV() const;
    int getQCRYS() const;
    int getDELTA4() const;
    int getFROZAT() const;
    int getNDEGF() const;
    int getNSTEP() const;
    int getNSAVC() const;
    int getNPRIV() const;
    const double* getPbc() const;
    const char* getTITLE() const;
    int getNTITLE() const;
    const int* getICNTRL() const;
    const char* getHDR() const;
    
    virtual ~DCD();
};

#endif // DCD_HPP_INCLUDED