/*
	File: libMC33.c
	Programmed by: David Vega
	February 2018
	updated by David Vega
	June 2019.
	February 2020
	August 2021
	September 2023
	February 2026
	March 2026
	Only this file must be included in a project to compile the library MC33
*/

/*****************************************************************************
You can change the following lines before compiling the library: */
#ifndef DEFAULT_SURFACE_COLOR
#define DEFAULT_SURFACE_COLOR 0xff5c5c5c /* RGBA 0xAABBGGRR: red 92 green 92 blue 92 (grey) */
#endif
#ifndef MC33_NORMAL_NEG
#define MC33_NORMAL_NEG 0 /* If it is 1, the front and back surfaces are exchanged. */
#endif
#ifndef USE_INTERNAL_SIGNBIT
#define USE_INTERNAL_SIGNBIT 1 /* definition of signbf function, see marching_cubes_33.c */
#endif
#ifndef USE_MM_RSQRT_SS
#define USE_MM_RSQRT_SS 1 /* definition of invSqrt function, see marching_cubes_33.c */
#endif
/*****************************************************************************/

#include "../include/marching_cubes_33.h"

#include "marching_cubes_33.c"
#include "MC33_util_grd.c"
