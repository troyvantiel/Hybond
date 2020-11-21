#include "stdafx.h"
#include "output.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;


//output of a bad volume file
//only gets called when the program tried to characterise more interactions than wnated which leads to very large volumes and runtimes.
void volOutput(string file){
		ofstream outputFile;
		outputFile.open(file + "p.csv");
		outputFile << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho" << endl;
		outputFile << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << endl;
		outputFile.close();

		outputFile.open(file + "m.csv");
		outputFile << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho" << endl;
		outputFile << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << endl;
		outputFile.close();

		outputFile.open(file + "t.csv");
		outputFile << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho difference" << endl;
		outputFile << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << "," << std::setw(20) << left << "na" << endl;
		outputFile.close();


}

//standard Bonder output method
void outputCube(double minx, double miny, double minz, double maxx, double maxy, double maxz, double res, string file, wfnData inputData,double cutoff,analysisBatch* batch,int makeCube)
{
	//std::cout << "Called Output cube" << std::endl;

	
	int dx = (maxx - minx) / res;
	int dy = (maxy - miny) / res;
	int dz = (maxz - minz) / res;

	int rdgLoops = (cutoff + 0.01) / 0.1 + 1;
	double *kinEngplus = new double[rdgLoops];
	double *potEngplus = new double[rdgLoops];
	double *totEngplus = new double[rdgLoops];
	double *pointsplus = new double[rdgLoops];
	double *elfplus = new double[rdgLoops];
	double *rhoplus = new double[rdgLoops];

	double *kinEngminus = new double[rdgLoops]; //set up of all the variable vectors for outputting to the various files
	double *potEngminus = new double[rdgLoops];
	double *totEngminus = new double[rdgLoops];
	double *pointsminus = new double[rdgLoops];
	double *elfminus = new double[rdgLoops];
	double *rhominus = new double[rdgLoops];

	//cout << makeCube << endl;

	for (size_t i = 0; i < rdgLoops; i++)
	{
		//std::cout << "just under the 1 print line" << std::endl;
		kinEngplus[i] =0;
		potEngplus[i] =0;
		totEngplus[i] = 0;
		pointsplus[i] = 0;
		elfplus[i] = 0;
		rhoplus[i] = 0;

		kinEngminus[i] = 0;
		potEngminus[i] = 0;
		totEngminus[i] = 0;
		pointsminus[i] = 0;
		elfminus[i] = 0;
		rhominus[i] = 0;
	}

	ofstream outputFile;
	if (makeCube)
		outputFile.open(file + ".cuberdg");
	ofstream outputFilerho;
	if (makeCube)
		outputFilerho.open(file + ".cuberhosign");

	//outputting to the cibe files for visulaising different parts of the molecule, such as the volume of the interaction that was looked at.
	if (makeCube)
	{
		outputFile << endl << endl;
		outputFile << "  " << inputData.nuc << "  " << minx << "  " << miny << "  " << minz << endl;
		outputFile << "  " << (dx-1) / makeCube + 1 << "  " << res * makeCube << "  0  0" << endl;
		outputFile << "  " << (dy-1) / makeCube + 1 << "  0  " << res * makeCube << "  0" << endl;
		outputFile << "  " << (dz-1) / makeCube + 1<< "  0  0  " << res * makeCube << endl;

		outputFilerho << endl << endl;
		outputFilerho << "  " << inputData.nuc << "  " << minx << "  " << miny << "  " << minz << endl;
		outputFilerho << "  " << (dx-1) / makeCube + 1 << "  " << res * makeCube << "  0  0" << endl;
		outputFilerho << "  " << (dy-1) / makeCube + 1 << "  0  " << res * makeCube << "  0" << endl;
		outputFilerho << "  " << (dz-1) / makeCube + 1 << "  0  0  " << res * makeCube << endl;
		//std::cout << "inputData.nuc: " << inputData.nuc << std::endl;
		for (size_t i = 0; i < inputData.nuc; i++)
		{
			outputFile << "  " << inputData.type[i] + 1 << "  " << inputData.type[i] + 1 << "  " << inputData.x[i] << "  " << inputData.y[i] << "  " << inputData.z[i] << endl;
			outputFilerho << "  " << inputData.type[i] + 1 << "  " << inputData.type[i] + 1 << "  " << inputData.x[i] << "  " << inputData.y[i] << "  " << inputData.z[i] << endl;
		}
	}





	double output,rho;
	int imod,jmod,kmod;
	for (size_t i = 0; i < dx; i++)
	{
		//std::cout << "Got further down" << std::endl;
		if(makeCube)
			imod = i % makeCube;
		//printf("grid is %f percent complete\n", ((i * 100.0) * (1.0/xs) ));
		for (size_t j = 0; j < dy; j++)
		{
			if (makeCube)
				jmod = j % makeCube;
			for (size_t k = 0; k < dz; k++)
			{
				if (makeCube)
					kmod = k % makeCube;

				output = (*batch).RDG_rho(minx + res * i, miny + res * j, minz + res * k, &rho);
				if (makeCube && !imod && !jmod && !kmod)
				{
					outputFile << "  " << output;
					outputFilerho << " " << rho;
					if ((k/makeCube) % 6 == 5)
					{
						outputFile << endl;
						outputFilerho << endl;
					}
				}

				if (output < cutoff)
				{
					double energy[3] ;
					(*batch).AKinEng(minx + res * i, miny + res * j, minz + res * k,energy);
					double elf = (*batch).elf();
					//double info = (*batch).Information(minx + res * i, miny + res * j, minz + res * k);
					//double gosh = (*batch).goshEntropy(minx + res * i, miny + res * j, minz + res * k);
					//double fisher = (*batch).elf(minx + res * i, miny + res * j, minz + res * k);
					for (size_t m = 0; (m < rdgLoops); m++)
					{
						if (output <= 0.1 * m)
						{
							if (rho > 0)
							{
								kinEngplus[m] += energy[0];
								potEngplus[m] += energy[1];
								totEngplus[m] += energy[2];
								//kinEngplus[m] += info;
								//potEngplus[m] += gosh;
								//totEngplus[m] += fisher;
								elfplus[m] += elf;
								pointsplus[m]++;
								rhoplus[m] += rho;
							}
							else
							{
								kinEngminus[m] += energy[0];
								potEngminus[m] += energy[1];
								totEngminus[m] += energy[2];
								//kinEngminus[m] += info;
								//potEngminus[m] += gosh;
								//totEngminus[m] += fisher;
								pointsminus[m]++;
								elfminus[m] += elf;
								rhominus[m] += rho;
							}


						}
					}
				}

			}
			if (makeCube && !imod && !jmod)
			{
				outputFile << endl;
				outputFilerho << endl;
			}

		}
	}
	for (size_t i = 0; i < rdgLoops; i++)
	{
		kinEngplus[i] *= res*res*res;
		potEngplus[i] *= res*res*res;
		totEngplus[i] *= res*res*res;
		elfplus[i] *= res*res*res;
		rhoplus[i] *=res*res*res;

		kinEngminus[i] *= res*res*res;
		potEngminus[i] *= res*res*res;
		totEngminus[i] *= res*res*res;
		elfminus[i] *= res*res*res;
		rhominus[i] *= res*res*res;
	}

	if (makeCube)
	{
		outputFile.close();
		outputFilerho.close();
	}
	//outputting the positive interactions to a file
	if (pointsplus[rdgLoops - 1] != 0)
	{
		outputFile.open(file + "p.csv");
		cout << "positive interaction" << endl;
		cout << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho" <<  endl;
		outputFile << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho" << endl;
		for (size_t i = 1; i < rdgLoops; i++)
		{
			cout << std::setw(20) << left << i * 0.1 << "," << std::setw(20) << left << pointsplus[i] << "," << std::setw(20) << left << pointsplus[i] * res * res * res << "," << std::setw(20) << left << kinEngplus[i] << "," << std::setw(20) << left << potEngplus[i] << "," << std::setw(20) << left << totEngplus[i] << "," << std::setw(20) << left << elfplus[i] << "," << std::setw(20) << left << rhoplus[i] << endl;

			outputFile << std::setw(20) << left << i * 0.1 << "," << std::setw(20) << left << pointsplus[i] << "," << std::setw(20) << left << pointsplus[i] * res * res * res << "," << std::setw(20) << left << kinEngplus[i] << "," << std::setw(20) << left << potEngplus[i] << "," << std::setw(20) << left << totEngplus[i] << "," << std::setw(20) << left << elfplus[i] << "," << std::setw(20) << left << rhoplus[i] << endl;
		}
		outputFile.close();
	}

	//outputting the negative(minus) interactions to a file
	if (pointsminus[rdgLoops - 1] != 0)
	{
		outputFile.open(file + "m.csv");
		cout << "negative interaction" << endl;
		cout << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho" << endl;

		outputFile << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho" << endl;
		for (size_t i = 1; i < rdgLoops; i++)
		{
			cout << std::setw(20) << left << i * 0.1 << "," << std::setw(20) << left << pointsminus[i] << "," << std::setw(20) << left << pointsminus[i] * res * res * res << "," << std::setw(20) << left << kinEngminus[i] << "," << std::setw(20) << left << potEngminus[i] << "," << std::setw(20) << left << totEngminus[i] << "," << std::setw(20) << left << elfminus[i] << "," << std::setw(20) << left << rhominus[i] << endl;

			outputFile << std::setw(20) << left << i * 0.1 << "," << std::setw(20) << left << pointsminus[i] << "," << std::setw(20) << left << pointsminus[i] * res * res * res << "," << std::setw(20) << left << kinEngminus[i] << "," << std::setw(20) << left << potEngminus[i] << "," << std::setw(20) << left << totEngminus[i] <<","<< std::setw(20) << left << elfminus[i] << "," << std::setw(20) << left << rhominus[i] << endl;
		}
		outputFile.close();
	}

	//outputting the total interactions to a file
	if (pointsplus[rdgLoops - 1] != 0 && pointsminus[rdgLoops - 1] != 0)
	{
		outputFile.open(file + "t.csv");
		cout << "total interaction" << endl;
		cout << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho difference" << endl;
		outputFile << std::setw(20) << left << "RDG" << std::setw(21) << left << ",number of points" << std::setw(21) << left << ",volume" << std::setw(21) << left << ",abr kinetic energy" << std::setw(21) << left << ",abr potental energy" << std::setw(21) << left << ",abr total energy" << std::setw(21) << left << ",elf" << std::setw(21) << left << ",rho difference" << endl;
		for (size_t i = 1; i < rdgLoops; i++)
		{
			cout << std::setw(20) << left << i * 0.1 << "," << std::setw(20) << left << pointsplus[i] + pointsminus[i] << "," << std::setw(20) << left << (pointsminus[i] + pointsplus[i]) * res * res * res << "," << std::setw(20) << left << kinEngplus[i] + kinEngminus[i] << "," << std::setw(20) << left << potEngplus[i] + potEngminus[i] << "," << std::setw(20) << left << totEngplus[i] + totEngminus[i] << "," << std::setw(20) << left << elfplus[i] + elfminus[i] << "," << std::setw(20) << left << rhoplus[i] + rhominus[i] << endl;
			outputFile << std::setw(20) << left << i * 0.1 << "," << std::setw(20) << left << pointsplus[i] + pointsminus[i] << "," << std::setw(20) << left << (pointsminus[i] + pointsplus[i]) * res * res * res << "," << std::setw(20) << left << kinEngplus[i] + kinEngminus[i] << "," << std::setw(20) << left << potEngplus[i] + potEngminus[i] << "," << std::setw(20) << left << totEngplus[i] + totEngminus[i] << "," << std::setw(20) << left << elfplus[i] + elfminus[i] << "," << std::setw(20) << left << rhoplus[i] + rhominus[i] << endl;
		}
		outputFile.close();
	}


	//releasing the vectors of data
	delete [] kinEngplus;
	delete [] potEngplus;
	delete [] totEngplus;
	delete [] pointsplus;
	delete [] rhoplus;
	delete [] elfplus;

	delete [] kinEngminus;
	delete [] potEngminus;
	delete [] totEngminus;
	delete [] pointsminus;
	delete [] elfminus;
	delete [] rhominus;
}
