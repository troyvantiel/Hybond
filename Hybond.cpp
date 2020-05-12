#include <Python.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pyhelp.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <istream>
#include <vector>



using namespace std;

int main()
{
	//Py_SetProgramName(argv[0]);
	//Py_SetPythonHome("C:/usr/lib/python3.6");

	//cout << "Test" << endl;
	Py_Initialize();
	cout << "Python initialized: " << Py_IsInitialized() << endl;
	PyRun_SimpleString("print ('Hello World')");

	Py_Finalize();

	
	char filename[] = "test.py";
	FILE* fp;

	Py_Initialize();
	fp = _Py_fopen(filename, "r");
	PyRun_SimpleFile(fp,filename);

	Py_Finalize();

	cout << "This is just before the helper class code" << endl;
	CPyInstance pyInstance;
	PyRun_SimpleString("print('Hello World from embedded python!!!')");

	return 0;
} 