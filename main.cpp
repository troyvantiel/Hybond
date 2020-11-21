#include "Compat.hpp"
#include "stdafx.h"
#include "analize.h"
#include "readwfn.h"
#include "iface.h"
#include "output.h"
#include <string>
#include <string.h>
#include <iostream>
#include <thread>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <fstream>
#include <argp.h>

const int SIZE = 100;

struct arguments
{
	char* type;
	char *inputFile="input.xyz";
	double res=0.02;     
	double cutoff=0.3;
	char *output="output"; 
	int cubesize=1;
	double x1=0,x2=1,y1=0,y2=1,z1=0,z2=1;
	int atom1=0,atom2=1,atom3=2,atom4=3;
	//char* dir=".";
	char* configFile="config";
};


static struct argp_option options[] =
{
	{"input",'i',"FILENAME",0,"The input xyz file"},
	{"resolution",'r',"VOXELSIZE",0,"The length in amstrongs of a single voxel in the intergration grid"},
	{"cutoff",'c',"RDGCUTTOFF",0,"The maximium value for RDG"},
	{"output",'o',"OUTPUTPREFIX",0,"The file prefix for all output files"},
	{"cubesize",'q',"CUBESCALE",0,"Reduces the resolution for output cubefiles by the given factor, 0 will prvent cube files from being written"},
	{"x1",1,"XSTART",0,"Either the x cordanate to start for point mode or the x cordante for one corner in cube mode"},
	{"y1",2,"YSTART",0,"Either the y cordanate to start for point mode or the y cordante for one corner in cube mode"},
	{"z1",3,"ZSTART",0,"Either the z cordanate to start for point mode or the z cordante for one corner in cube mode"},
	{"x2",4,"XEND",0,"The x cordanate for the second corner in cube mode"},
	{"y2",4,"YEND",0,"The y cordanate for the second corner in cube mode"},
	{"z2",4,"ZEND",0,"The z cordanate for the second corner in cube mode"},
	{"atom1",'1',"ATOMNUMBER",0,"The first atom for use in line, trinagle or quad mode"},
	{"atom2",'2',"ATOMNUMBER",0,"The first atom for use in line, trinagle or quad mode"},
	{"atom3",'3',"ATOMNUMBER",0,"The first atom for use in  trinagle or quad mode"},
	{"atom4",'4',"ATOMNUMBER",0,"The first atom for use in quad mode"},
	//{"directory",'d',"FOLDER",0,"The folder to put output files in"},
	{"config",'f',"CONFIGFILE",0,"A file containing all of the options to run bonder, only in input mode"},
	{0}
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	arguments *argument = (arguments*)(state->input);

	switch (key)
	{
		case 'i':
			argument->inputFile = arg;
			break;
		case 'r':
			argument->res = std::stod(arg);
			break;
		case 'c':
			argument->cutoff = std::stod(arg);
			break;
		case 'o':
			argument->output = arg;
			break;
		case 'q':
			argument->cubesize = std::stoi(arg);
                        break;
		case 1:
			argument->x1 = std::stod(arg);
                        break;
		case 2:
                        argument->y1 = std::stod(arg);
                        break;
		case 3:
                        argument->z1 = std::stod(arg);
                        break;
		case 4:
                        argument->x2 = std::stod(arg);
                        break;
		case 5:
                        argument->y2 = std::stod(arg);
                        break;
		case 6:
                        argument->y1 = std::stod(arg);
                        break;
		case '1':
                        argument->atom1 = std::stoi(arg);
                        break;
		case '2':
                        argument->atom2 = std::stoi(arg);
                        break;
		case '3':
                        argument->atom3 = std::stoi(arg);
                        break;
		case '4':
                        argument->atom4 = std::stoi(arg);
                        break;
		//case 'd':
                //        argument->dir = arg;
                //        break;
		case 'f':
			argument->configFile=arg;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 1)
			{
				argp_usage(state);
			}
			else
			{
				argument->type = arg;
			}
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 1)
			{
				argp_usage (state);
			}
			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static char args_doc[] = "type";
static char doc[] = "Bonder, a program designed to identify and map non covalent interactions";
static struct argp argp = {options, parse_opt, args_doc, doc};


wfnData* init(std::string file)
{
	wfnData *inputFile = readFile(file);
	return inputFile;
}

void drawline(int a, int b, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube, int outputfilemod)
{
	//std::cout << "int a is :" << a << " and int b is: " << b << std::endl;
	double lowX = (*inputFile).x[a];
	double lowY = (*inputFile).y[a];
	double lowZ = (*inputFile).z[a];

	double highX = (*inputFile).x[b];
	double highY = (*inputFile).y[b];
	double highZ = (*inputFile).z[b];


	//std::cout << "testing line with atom coords of:" <<std::endl;
	//std::cout << "LowX: " << lowX << " LowY: " << lowY << " LowZ: " << lowZ << std::endl;         //these comments are added in by me
	//std::cout << "HighX: " << highX << " HighY: " << highY << " HighZ: " << highZ << std::endl;


	//printf("testing line between %d and %d\n",a,b);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	analysisBatch* batch = new analysisBatch(*inputFile);
	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	bool sucsess = false;

	//printf("%d\n",reps);
	int flips = 1;
	//printf("size of reps which analizes a point everytime is %d\n", reps);

	double lowest = 9999;

	for (size_t i = 0; i < reps; i++)
	{
		int k = reps / 2 + flips * i / 2;
		flips *= -1;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured < lowest)
			lowest = mesured;
		if (mesured <= cutoff)
		{
			//std::cout << "found point" << std::endl;

			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);
			outputfile += "Frame:" + std::to_string(outputfilemod);
			analize.anilizePoint(0, 0, 0, 0, size, size, cutoff, &sucsess, inputFile, outputfile, batch,makeCube);
			
			break;
		}
	}
	//std::cout << lowest << std::endl;

	delete batch;
}

void drawtrig(int a, int b,int c, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{

	analysisBatch* batch = new analysisBatch(*inputFile);
	double highX = ((*batch).atomx(a) + (*batch).atomx(b))/2;
	double highY = ((*batch).atomy(a) + (*batch).atomy(b))/2;
	double highZ = ((*batch).atomz(a) + (*batch).atomz(b))/2;

	double lowX = (*batch).atomx(c);
	double lowY = (*batch).atomy(c);
	double lowZ = (*batch).atomz(c);

	printf("testing line between centre of %d and %d and %d\n",a,b,c);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	bool sucsess = false;

	printf("%d\n",reps);
	for (size_t i = 0; i < reps; i++)
	{
		int k = i;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured <= cutoff)
		{

			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);

			analize.anilizePoint(0, 0, 0, 0, size, size, cutoff, &sucsess, inputFile, outputfile, batch,makeCube);
			break;
		}
	}


	delete batch;
}

void drawquad(int a, int b,int c,int d, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{

	analysisBatch* batch = new analysisBatch(*inputFile);
	double highX = ((*batch).atomx(a) + (*batch).atomx(b))/2;
	double highY = ((*batch).atomy(a) + (*batch).atomy(b))/2;
	double highZ = ((*batch).atomz(a) + (*batch).atomz(b))/2;

	double lowX = ((*batch).atomx(c)+(*batch).atomx(d))/2;
	double lowY = ((*batch).atomy(c)+(*batch).atomx(d))/2;
	double lowZ = ((*batch).atomz(c)+(*batch).atomx(d))/2;

	printf("testing line between centre of %d and %d and %d\n",a,b,c);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	bool sucsess = false;

	printf("%d\n",reps);
	for (size_t i = 0; i < reps; i++)
	{
		int k = i;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured <= cutoff)
		{

			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);

			analize.anilizePoint(0, 0, 0, 0, size, size, cutoff, &sucsess, inputFile, outputfile, batch,makeCube);
			break;
		}
	}


	delete batch;
}
struct pdrawArgs
{
	int a; int b; double res; double cutoff; std::string outputfile; int size; wfnData* inputFile; int makeCube;
	pdrawArgs(int A, int B, double Res, double Cutoff, std::string Outputfile, int Size, wfnData* InputFile,int MakeCube)
	{
		makeCube = MakeCube;
		a = A;
		b = B;
		res = Res;
		cutoff = Cutoff;
		outputfile = Outputfile;
		size = Size;
		inputFile = InputFile;
	}
};

void pDrawline(void *input)
{
	pdrawArgs* data = (pdrawArgs*)input;
	drawline((*data).a, (*data).b, (*data).res, (*data).cutoff, (*data).outputfile, (*data).size, (*data).inputFile, (*data).makeCube, 0);
	//pthread_exit(NULL);
	delete data;
}

void runAll(double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube)
{



	//set up threadpool
	boost::asio::io_service ioService;
	boost::thread_group threadpool;
	std::unique_ptr<boost::asio::io_service::work> work(new boost::asio::io_service::work(ioService));

	int numOfThreads = std::thread::hardware_concurrency();
	std::cout << numOfThreads << std::endl;
	if (numOfThreads > 24)
		numOfThreads = 24;
	for(int i = 0; i< numOfThreads; i++)
		threadpool.create_thread(boost::bind(&boost::asio::io_service::run, &ioService));


	for (size_t i = 0; i < (*inputFile).nuc; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			//997 is a large prime, will break on protens of more than 2 million atoms 
			pdrawArgs *lineData;
			lineData = new pdrawArgs(i, j, res, cutoff, outputfile, size, inputFile, makeCube);
			ioService.post(boost::bind(pDrawline, (void *)lineData));

		}
	}

	work.reset();
	threadpool.join_all();
	ioService.stop();

}

std::vector<std::string> readFileLines(const char* filename)
{
	std::vector<std::string> file;
	std::ifstream input(filename);
	std::string line;
	int i = 0;
	while (getline(input, line)){
		file.push_back(line);
	}
	return file;
}
void useInputFile(char* filename)
{

	std::fstream inputFileTest(filename);
	if(!inputFileTest)
	{
		std::cout << "input file not found" << std::endl;
		return;
	}
	std::vector<std::string> lines;
	lines = readFileLines(filename);
	int lineNum = lines.size();
	if(lineNum == 0)
	{
		std::cout << "the input file needs text" <<  std::endl;
		//printf("bonder h for help\n");
		return;
	}

	if(lineNum == 1)
	{
		std::cout << "please select option in the input file and the name of the wfn file";
		return;
	}

	wfnData *inputFile = 0;
	try
	{
		inputFile = init(lines[1]);

	}
	catch (const std::invalid_argument& ia) 
	{
		std::cout << "error in parssing XYZ data" << std::endl;
		return;
	}

	if(lines[0] == "p")
	{
		if (lineNum != 9)
		{
			std::cout << "error in parsing input file\npoint file format is:\np\nwfn file\nx\ny\nz\nrdg cutoff\nres\noutput fiel name\noutput cube file"<< std::endl;
			return;
		}
		bool sucsess;
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		try
		{
			analize.setUpAnalysisBatch( std::stod(lines[2]), std::stod(lines[3]), std::stod(lines[4]), std::stod(lines[5]),batch);
			printf("%f \n", (*batch).RDG(std::stod(lines[2]), std::stod(lines[3]), std::stod(lines[4])));
			analize.anilizePoint(0, 0, 0, 0, SIZE, SIZE, std::stod(lines[6]), &sucsess, inputFile, lines[7],batch, std::stoi(lines[8]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}

		if (sucsess)
		{
			printf("point given is in region\n");
		}
		else
		{
			printf("point given is not in region\n");
		}
		return;

	}

	if (lines[0] == "l")
	{
		if (lineNum != 8)
		{
			std::cout << "error in parsing input file\nline file format is:\nl\nwfn file\natom1\natom2\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}   
		try
		{
			drawline(std::stoi(lines[2]), std::stoi(lines[3]), std::stod(lines[4]), std::stod(lines[5]), lines[6], SIZE, inputFile, std::stoi(lines[7]), 0);
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;

	}
	//letter file 1 2 res cutoff
	if (lines[0] == "t")
	{
		if (lineNum != 9)
		{
			std::cout << "error in parsing input file\ntriangle file format is:\nt\nwfn file\natom1\natom2\natom3\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}

		try
		{
			drawtrig(std::stoi(lines[2]), std::stoi(lines[3]),std::stoi(lines[4]), std::stod(lines[5]), std::stod(lines[6]), lines[7], SIZE, inputFile, std::stoi(lines[8]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;

	}

	if (lines[0] == "q")
	{
		if (lineNum != 10)
		{
			std::cout << "error in parsing input file\ntriangle file format is:\nt\nwfn file\natom1\natom2\natom3\natom4\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}

		try
		{
			drawquad(std::stoi(lines[2]), std::stoi(lines[3]),std::stoi(lines[4]),std::stoi(lines[5]), std::stod(lines[6]), std::stod(lines[7]), lines[8], SIZE, inputFile, std::stoi(lines[9]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;

	}
	//letter file res cutoff output
	if (lines[0] == "a")
	{
		if (lineNum != 6)
		{
			std::cout << "error in parsing input file\nall bonds file format is:\na\nwfn file\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}


		try
		{
			runAll(std::stod(lines[2]), std::stod(lines[3]), lines[4], SIZE, inputFile, std::stoi(lines[5]));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;
	}

	//letter file minx miny minz maxx maxy maxz res outputFile
	if (lines[0] == "g")
	{
		if (lineNum != 10)
		{
			std::cout << "error in parsing input file\ngrid file format is:\ng\nwfn file\nlow x\nlow y\n low z\n high x\n high y \n high z\nres\noutput file name" << std::endl;
			return;
		}

		analysisBatch* batch = new analysisBatch(*inputFile);
		try
		{

			outputCube(std::stod(lines[2]), std::stod(lines[3]), std::stod(lines[4]), std::stod(lines[5]), std::stod(lines[6]), std::stod(lines[7]), std::stod(lines[8]), lines[9], *inputFile, 1.0, batch, 1);
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		printf("done");
		return;
	}

}


void bond(int argc, char* argv[], std::string newfile, int outputfilemod) //was int main
{
	struct arguments arguments;
	argp_parse (&argp, argc, argv, 0, 0, &arguments);


	if (arguments.type[0] == 'f')
	{
		useInputFile(arguments.configFile);
		//return 0;
	}

	wfnData *inputFile = 0;
		try
		{

			inputFile = init(newfile);

			//inputFile = init(arguments.inputFile);

		}
		catch (const std::invalid_argument& ia) 
		{
			std::cout << "error in parssing molcule data" << std::endl;
			//return 1;
		}

	//std::cout << "data read" << std::endl;
	//letter file x y z res cutoff
	if (arguments.type[0] == 'p')
	{
		bool sucsess;
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		analize.setUpAnalysisBatch( arguments.x1 , arguments.y1, arguments.z1, arguments.res,batch);
		analize.anilizePoint(0, 0, 0, 0, SIZE, SIZE, arguments.cutoff, &sucsess, inputFile, arguments.output, batch, arguments.cubesize);

		if (sucsess)
		{
			printf("point given is in region\n");
		}
		else
		{
			printf("point given is not in region\n");
		}
		//return 0;
	}

	//letter file 1 2 res cutoff
	if (arguments.type[0] == 'l')
	{
		std::cout << "Line Calculation Selected" << std::endl;
		drawline(arguments.atom1, arguments.atom2, arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize, outputfilemod);
		//return 0;

	}

	//letter file 1 2 res cutoff
	if (arguments.type[0] == 't')
	{
				drawtrig(arguments.atom1, arguments.atom2,arguments.atom3, arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		//return 0;

	}

	if (arguments.type[0] == 'q')
	{
		drawquad(arguments.atom1, arguments.atom2,arguments.atom3,arguments.atom4, arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		//return 0;
	}


	//letter file res cutoff output
	if (arguments.type[0] == 'a')
	{
		runAll(arguments.res, arguments.cutoff, arguments.output, SIZE, inputFile, arguments.cubesize);
		//return 0;
	}

	//letter file minx miny minz maxx maxy maxz res outputFile
	if (arguments.type[0] == 'g')
	{
		analysisBatch* batch = new analysisBatch(*inputFile);
		outputCube(arguments.x1, arguments.y1, arguments.z1, arguments.x2, arguments.y2, arguments.z2, arguments.res, arguments.output, *inputFile, 1.0, batch, arguments.cubesize);
		printf("done");
		//return 0;
	}


	//bool sucsess;
	//anilizePoint(0, 0, 0, 0, 2000, 2000, 2501, &sucsess);
	//printf("%d", sucsess);
	//printf("bonder h for help\n");
	//return 0;
}

