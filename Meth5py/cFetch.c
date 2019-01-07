#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "methStruct.h"

//http://www.codingunit.com/c-tutorial-binary-file-io

void fetch(char* fileName, long seekStart, int blocks, 
	int minCov, int maxCov, unsigned short *c, unsigned short *ct, unsigned char *con) {
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
	int i;
	double tmp;
	for(i=0; i<blocks; i++) {
		if(structArray[i].ct == 65535) {
			c[i] = 0;
			ct[i] = 0;
		} else if(structArray[i].ct < (unsigned short)minCov) {
			c[i] = 0;
			ct[i] = 0;
		} else if(structArray[i].ct > (unsigned short)maxCov) {
			c[i] = 0;
			ct[i] = 0;
		} else {
			c[i] = structArray[i].c;
			ct[i] = structArray[i].ct;
		}
		con[i] = (unsigned char) structArray[i].con;
	}
	free(structArray);
}

static PyObject* cFetch(PyObject* self, PyObject* args) {
	// Pass back a tuple of arrays
	char* fileName = NULL;
	long seekStart = 0;
	int blocks = 0;
	int minCov = 0;
	int maxCov = 0;
	if (!PyArg_ParseTuple(args, "sliii", &fileName, &seekStart, &blocks,
		&minCov, &maxCov)) {
		return NULL;
	}
	npy_intp dim[1] = {blocks};
	PyObject* outC = PyArray_SimpleNew(1, dim, NPY_UINT16);
	PyObject* outCT = PyArray_SimpleNew(1, dim, NPY_UINT16);
	PyObject* outContext = PyArray_SimpleNew(1, dim, NPY_UINT8);
	unsigned short *c = (unsigned short*)PyArray_DATA(outC);
	unsigned short *ct = (unsigned short*)PyArray_DATA(outCT);
	unsigned char *con = (unsigned char*)PyArray_DATA(outContext);
	//PyObject* outMeth = PyList_New((Py_ssize_t) blocks);
	//PyObject* outContext = PyList_New((Py_ssize_t) blocks);
	fetch(fileName, seekStart, blocks, minCov, maxCov, c, ct, con);
	return Py_BuildValue("(NNN)", outC, outCT, outContext);
}

static PyMethodDef FetchMethods[] = {
	{"cFetch", cFetch, METH_VARARGS, "Implements fetch in C"},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcFetch(void) {
	(void) Py_InitModule("cFetch", FetchMethods);
	import_array();
}
