#include <cstdlib>
#include <cstring>
#include <signal.h>
#include <stdlib.h>

#include <fstream>
#include <stdio.h>
//<>
#include <iostream>

#include "dcd_r.hpp"

using namespace std;

DCD_R::DCD_R(const char filename[])
{
	dcdf.exceptions(std::ifstream::failbit);
	try
	{
		//dcdf
		dcdf.open(filename, ios::in|ios::binary);
	}
	catch(std::ifstream::failure e)
	{
		cerr << "Exception opening/reading file '" << filename << "' : " << endl;
		cerr << "Please check the path of the file and if it exists." << endl;
	}

	dcd_first_read = true;
} 

void DCD_R::alloc()
{
	X = new float[NATOM];
	Y = new float[NATOM];
	Z = new float[NATOM];
	pbc[0] = pbc[1] = pbc[2] = pbc[3] = pbc[4]= pbc[5] = 0.0;
}
void DCD_R::read_header()
{
	unsigned int fortcheck1, fortcheck2;

	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));    //FORT CHECK NEEDS UNCOMMENT
	dcdf.read((char*)HDR,sizeof(char)*4);
	dcdf.read((char*)ICNTRL,sizeof(int)*20);
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));			//FORT CHECK NEEDS UNCOMMENT
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	HDR[4] = '\0';
	NFILE = ICNTRL[0];
	NPRIV = ICNTRL[1];
	NSAVC = ICNTRL[2];	
	NSTEP = ICNTRL[3];
	NDEGF = ICNTRL[7];
	FROZAT = ICNTRL[8];
	DELTA4 = ICNTRL[9];
	QCRYS = ICNTRL[10];
	CHARMV = ICNTRL[19];

	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));    //FORT CHECK NEEDS UNCOMMENT
	dcdf.read((char*)&NTITLE, sizeof(int));

	if(NTITLE == 0)
	{
		TITLE = new char[80+1];
		TITLE[0] = '\0';
	}
	else
	{
		TITLE = new char[NTITLE*80+1];
		for(int i = 0; i < NTITLE; i++)
		{
			dcdf.read((char*)&TITLE[i*80], sizeof(char)*80);
		}
		TITLE[NTITLE*80] = '\0';
	}
	dcdf.read((char*)&fortcheck2, sizeof(unsigned int));            //FORT CHECK NEEDS UNCOMMENT
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	//Reading number of atoms
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)&NATOM,sizeof(int));
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2); //FORT CHECK UNCOMMENT

	LNFREAT = NATOM - FROZAT;
	if(LNFREAT != NATOM)
	{
		FREEAT = new int[LNFREAT];
		dcdf.read((char*)&fortcheck1, sizeof(unsigned int));
		dcdf.read((char*)FREEAT, sizeof(int)*LNFREAT);
		dcdf.read((char*)&fortcheck2, sizeof(unsigned int));
		checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2); //FORT CHECK UNCOMMENT
	}
	alloc();
}

void DCD_R::read_oneFrame()
{
	//try
	//{
		unsigned int fortcheck1, fortcheck2 =0;

	int siz = (dcd_first_read) ? NATOM : LNFREAT;
	//int siz = LNFREAT;
	float *tmpX = new float[siz];
	float *tmpY = new float[siz];
	float *tmpZ = new float[siz];

	if(QCRYS)
	{
		//cout<<"QCHRYS count was greater than 0" << endl;
		dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
		dcdf.read((char*)pbc,sizeof(double)* 6);
		dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
		checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2); //FORTCHECK UNCOMMENT
	}
		//X
	//cout <<"Reading co-ordinates" << endl;
	//cout <<"***fortcheck1***"<< endl;
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	cout<< "fortcheck1: " << fortcheck1<<endl;
	//cout << "*****reading x*****" <<endl;
	dcdf.read((char*)tmpX,sizeof(float)*siz);	
	//cout << "***fortcheck2***"<<endl;
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
	//cout<< "fortcheck2: " << fortcheck2 << endl;
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);
	//cout<<"finished getting x co-ords"<<endl;
	//Y
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)tmpY,sizeof(float)*siz);	
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	//Z
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)tmpZ,sizeof(float)*siz);	
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);
	//cout << "finished reading co-ordinates" <<endl;
	


	if(dcd_first_read)
	{
		//cout << "memcpy on first read" << endl;
		memcpy(X, tmpX, NATOM*sizeof(float));
		memcpy(Y, tmpY, NATOM*sizeof(float));
		memcpy(Z, tmpZ ,NATOM*sizeof(float));
	}
	else
	{
		if(LNFREAT != NATOM)
		{
			//cout << "before the for loop to reading xyz" <<endl;
			for(int i =0; i < siz; i++)
			{
				X[FREEAT[i]-1] = tmpX[i];
				Y[FREEAT[i]-1] = tmpY[i];
				Z[FREEAT[i]-1] = tmpZ[i];
			}
		}
		else
		{
			//cout<<"memcpy else condition"<<endl;
			memcpy(X, tmpX, NATOM*sizeof(float));
			memcpy(Y, tmpY, NATOM*sizeof(float));
			memcpy(Z, tmpZ ,NATOM*sizeof(float));
		}
	}

	if(dcd_first_read)
		dcd_first_read=false;
	//cout << "delete temp arrays" <<endl;
	delete[] tmpX;
	delete[] tmpY;
	delete[] tmpZ;
	}
	//catch(const std ::exception& e)
	//{
	//	cout << e.what() <<endl;
	//}
	
//}
void DCD_R::printHeader() const
{
    int i;
    
    cout << "HDR :\t" << HDR << endl;
    

    cout << "ICNTRL at index 0 = Number of Frames" << endl;
    cout << "ICNTRL at index 1 = if restart, total number of frames before first print " << endl;
    cout << "ICNTRL at index 2 = frequency of writing dcd" << endl;
    cout << "ICNTRL at index 3 = Number of steps ; note: steps/freq = no of frames" << endl;
    cout << "ICNTRL at index 7 = Number of degrees of freedom" << endl;
    cout << "ICNTRL at index 8 = Is number of free atoms" << endl;
    cout << "ICNTRL at index 9 = timestep in AKMA units" << endl;
    cout << "ICNTRL at index 10 = is 1 if CRYSTAL used" << endl;
    cout << "ICNTRL at index 19 = CHARMM version" << endl;

    cout << "ICNTRL :\t";
    for(i=0;i<20;i++)
        cout << ICNTRL[i] << "\t" ;
    cout << endl;
    
    cout << "NTITLE :\t" << NTITLE << endl;
    cout << "TITLE :\t" << TITLE << endl;
    
    cout << "NATOM :\t" << NATOM << endl;
    cout << "LNFREAT :\t" << LNFREAT << endl;
    
}

DCD_R::~DCD_R()
{
    dcdf.close();
    
    delete[] TITLE;
    
    if (LNFREAT != NATOM)
        delete[] FREEAT;
    
    delete[] X;
    delete[] Y;
    delete[] Z;
}
