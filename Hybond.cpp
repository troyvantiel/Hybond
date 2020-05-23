//#define PY_SSIZE_T_CLEAN

//#include <Python.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "pyhelp.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <istream>
#include <vector>



using namespace std;

int main(int argc, char* argv[])
{
	/*
	int len = sizeof(argc)/1;
	cout << "Length of input array: " << len << endl;
	if(len > 4)
	{
		cout<< "Usage: DCD File input , Output file" << endl;
		return 1;
	}*/
		streampos size;
		char * memblock;

		ifstream fileBuffer(argv[0], ios::in|ios::binary|ios::ate);
		cout<< "Made ifstream buffer" << endl;
		ofstream outputBuffer("output.txt", ios::out|ios::binary);
		cout<< "Made ofstream buffer" << endl;

		char input[1024];
		char output[1024];

		if(fileBuffer.is_open())
		{
			size = fileBuffer.tellg();
			memblock = new char[size];
			fileBuffer.seekg(0, ios::beg);
			fileBuffer.read(memblock,size);
		}

		outputBuffer.write((char*)memblock, size);
		outputBuffer.close();
		fileBuffer.close();
		return 0;

}






//Following code was written to run python scripts to read dcd files but the methods couldnt
//be run through the C++ code. This means the DCD file could not be passed to the python reader
//and an output back to the C++ code could not be obtained.
/*
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
	

	cout << "Start of the get integer code" << endl;

	//CPyInstance pyInstance;

	CPyObject pName = PyUnicode_FromString("test");
	CPyObject pModu = PyImport_Import(pName);

	if(pModu)
	{
		CPyObject pFunc = PyObject_GetAttrString(pModu, "getInteger");
		if(PyCallable_Check(pFunc)==1)
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

int main(int argc , char *argv[])
{
	PyObject *pName, *pModule, *pFunc;
	PyObject *pArgs, *pValue;
	int i;

	if(argv < 3)
	{
		printf(stderr, "Usage: Call python file func name [args]\n");
		return 1;
	}

	Py_Initialize();
	pName = PyString_FromString(argv[1]);

	pModule = PyImport_Import(pName);
	Py_DECREF(pName);

	if(pModule != NULL) 
	{
		pFunc = PyObject_GetAttrString(pModule, argv[2]);

		if(pFunc && PyCallable_Check(pFunc))
		{
			pArgs = PyTuple_New(argc - 3);
			for(i = 0; i < argc - 3; ++i)
			{
				pValue = PyInt_FromLong(atoi(argv[i+3]));
				if(!pValue)
				{
					Py_DECREF(pArgs);
					Py_DECREF(pModule);
					printf(stderr, "Cannot convert argument\n");
					return 1;
				}
				PyTuple_SetItem(pArgs,i,pValue);
			}
			pValue = PyObject_CallObject(pFunc,pArgs);
			Py_DECREF(pArgs);
			if(pValue !=NULL)
			{
				printf("Result of call: %ld\n", PyInt_AsLong(pValue));
				Py_DECREF(pValue);
			}
			else
			{
				Py_DECREF(pFunc);
				Py_DECREF(pModule);
				PyErr Print();
				fprintf(stderr, "Call Failed\n");
				return 1;
			}
		}
		else
		{
			if(PyErr Occurred())
				PyErr_Print();
			fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
		}
		Py_XDECREF(pFunc);
		Py_DECREF(pModule);
	}
	else
	{
		PyErr_Print();
		fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
		return 1;
	}
	Py_Finalize();
	return 0;
}*/

