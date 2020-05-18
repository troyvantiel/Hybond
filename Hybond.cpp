#define PY_SSIZE_T_CLEAN

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

	/*
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
	*/

	cout << "Start of the get integer code" << endl;

	CPyInstance hInstance;

	CPyObject pName = PyUnicode_FromString("test");
	CPyObject pModule = PyImport_Import(pName);

	if(pModule)
	{
		CPyObject pFunc = PyObject_GetAttrString(pModule, "getInteger");
		if(pFunc && PyCallable_Check(pFunc))
		{
			CPyObject pValue = PyObject_CallObject(pFunc, NULL);

			printf("C: getInteger() = %ld\n", PyLong_AsLong(pValue));
		}
		else
		{
			printf("ERROR: function getInteger()\n");
		}

	}
	else
	{
		printf("ERROR: Module not imported\n");
	}
	return 0;
} 