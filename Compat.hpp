
//definitions of the methods go here but only the signatures
#include <string>
#include <string.h>
#include <iostream>
#include <thread>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <fstream>
#include <argp.h>
#include <readwfn.h>

static error_t parse_opt(int key, char *arg, struct argp_state *state);

void drawline(int a, int b, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube, int outputfilemod);

void drawtrig(int a, int b,int c, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube);

void drawquad(int a, int b,int c,int d, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube);

void pDrawline(void *input);

void runAll(double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,int makeCube);

std::vector<std::string> readFileLines(const char* filename);

void useInputFile(char* filename);

void bond(int argc, char *argv[], std::string newfile, int outputfilemod);

wfnData* init(std::string file);



