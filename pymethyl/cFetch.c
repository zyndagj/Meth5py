#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include "methStruct.h"

//http://www.codingunit.com/c-tutorial-binary-file-io

double* fetch(char* fileName, int seekStart, int blocks, 
	int minCov, int maxCov) {
	FILE* f=fopen(fileName,"rb");
	int fseekError = fseek(f, (long int) seekStart, SEEK_SET);
	// Check for errors in fseek
	if(fseekError) {
		printf("Seek error\n");
	}
	methStruct* structArray = malloc(sizeof(methStruct)*blocks);
	size_t readRet = fread(structArray, sizeof(methStruct), 
		(size_t) blocks, f);
	if(readRet != blocks) {
		printf("Read failed.\n");
	}
	fclose(f);
	double* outArray = malloc(sizeof(double)*blocks);
	int i;
	for(i=0; i<blocks; i++) {
		if(structArray[i].x == (unsigned short) 65535 && structArray[i].y == (unsigned short) 65535) {
			outArray[i] = -1;
		} else if(structArray[i].y < minCov) {
			outArray[i] = -1;
		} else if(structArray[i].y > maxCov) {
			outArray[i] = -1;
		} else {
			outArray[i] = (double) structArray[i].x / (double) structArray[i].y;
		}
	}
	free(structArray);
	return outArray;
}

static PyObject* cFetch(PyObject* self, PyObject* args) {
	char* fileName = NULL;
	int seekStart = 0;
	int blocks = 0;
	int minCov = 0;
	int maxCov = 0;
	if (!PyArg_ParseTuple(args, "siiii", &fileName, &seekStart, &blocks,
		&minCov, &maxCov)) {
		return NULL;
	}
	
	double* outArray = fetch(fileName, seekStart, blocks, minCov, maxCov);
	PyObject* outList = PyList_New((Py_ssize_t) blocks);
	int i;
	for(i=0; i<blocks; i++) {
		PyList_SetItem(outList, (Py_ssize_t) i, Py_BuildValue("d",outArray[i]));
	}
	free(outArray);
	return outList;
}

static PyMethodDef FetchMethods[] = {
	{"cFetch", cFetch, METH_VARARGS, "Implements fetch in C"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcFetch(void) {
	(void) Py_InitModule("cFetch", FetchMethods);
}
