#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include "dcd_r.hpp"

using namespace std;

DCD_R::DCD_R(const char filename[])
{
	dcdf.exceptions(std::ifstream::failbit);
	try
	{
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

	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)HDR,sizeof(char)*4);
	dcdf.read((char*)ICNTRL,sizeof(int)*20);
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
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

	dcdf.read((char*)&fortcheck1,sizeof(unisgned int));
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
	dcdf.read((char*)&fortcheck2, sizeof(unisgned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	//Reading number of atoms
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)&NATOM,sizeof(int));
	dcdf.read((char*)&fortcheck2,sizeof(unsigned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	LNFREAT = NATOM - FROZAT;
	if(LNFREAT != NATOM)
	{
		FREAT = new int[LNFREAT];
		dcdf.read((char*)&fortcheck1, sizeof(unsigned int));
		dcdf.read((char*)FREAT, sizeof(int)*LNFREAT);
		dcdf.read((char*)&fortcheck2, sizeof(unsigned int));
		checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);
	}
	alloc();
}

void DCD_R::read_oneFrame()
{
	unisgned int fortcheck1, fortcheck2;

	int siz = (dcd_first_read) ? NATOM : LNFREAT;

	float *tmpX = new float[siz];
	float *tmpY = new float[siz];
	float *tmpZ = new float[siz];

	if(QCRYS)
	{
		dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
		dcdf.read((char*)pbc,sizeof(double)*6);
		dcdf.read((char*)&fortcheck2,sizeof(unisgned int));
		checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);
	}

	//X
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)tmpX,sizeof(float)*siz);	
	dcdf.read((char*)&fortcheck2,sizeof(unisgned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	//Y
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)tmpY,sizeof(float)*siz);	
	dcdf.read((char*)&fortcheck2,sizeof(unisgned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	//Z
	dcdf.read((char*)&fortcheck1,sizeof(unsigned int));
	dcdf.read((char*)tmpZ,sizeof(float)*siz);	
	dcdf.read((char*)&fortcheck2,sizeof(unisgned int));
	checkFortranIOerror(__FILE__,__LINE__,fortcheck1,fortcheck2);

	if(dcd_first_read)
	{
		memcpy(X, tmpX, NATOM*sizeof(float));
		memcpy(Y, tmpY, NATOM*sizeof(float));
		memcpy(Z, tmpZ ,NATOM*sizeof(float));
	}
	else
	{
		if(LNFREAT != NATOM)
		{
			for(int i =0; i < siz; i++)
			{
				X[FREEAT[i]-1] = tmpX[i];
				Y[FREEAT[i]-1] = tmpY[i];
				Z[FREEAT[i]-1] = tmpZ[i];
			}
		}
		else
		{
			memcpy(X, tmpX, NATOM*sizeof(float));
			memcpy(Y, tmpY, NATOM*sizeof(float));
			memcpy(Z, tmpZ ,NATOM*sizeof(float));
		}
	}

	if(dcd_first_read)
		dcd_first_read=false;

	delete[] tmpX;
	delete[] tmpY;
	delete[] tmpZ;
}
void DCD_R::printHeader() const
{
    int i;
    
    cout << "HDR :\t" << HDR << endl;
    
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
