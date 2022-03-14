/*
	File: libMC33.c
	Programmed by: David Vega
	February 2018
	updated by David Vega
	June 2019.
	February 2020
	August 2021
	Only this file must be included in a project to compile the library MC33
*/

#define compiling_libMC33
/********************************CUSTOMIZING**********************************/
//The following line can be only changed before compiling the library:
//#define MC_Normal_neg // the front and back surfaces are exchanged.
/*****************************************************************************/

#include "../include/marching_cubes_33.h"

#ifndef GRD_orthogonal
extern void (*mult_Abf)(double (*)[3], float *, float *, int);
extern int DefaultColorMC;
#endif

#include "marching_cubes_33.c"
