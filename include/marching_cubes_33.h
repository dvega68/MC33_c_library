/*
	File: marching_cubes_33.h
	Programmed by: David Vega - dvega@uc.edu.ve
	           Javier Abache - jabache@uc.edu.ve
	March 2012
	updated by David Vega:
	December 2014,
	June 2015,
	November 2017.
	February 2018.
	June 2019.
	January 2020.
	February 2020
	July 2020
	February 2021
	August 2021
	December 2021
	January 2022
	September 2023
	This library is a modified version of the library described in:
	Vega, D., Abache, J., Coll, D., A Fast and Memory Saving Marching Cubes 33
	implementation with the correct interior test, Journal of Computer Graphics
	Techniques (JCGT), vol. 8, no. 3, 1–18, 2019
*/

#ifndef marching_cubes_33_h
#define marching_cubes_33_h
#if defined(compiling_libMC33)
#undef integer_GRD
#undef size_type_GRD
#undef GRD_orthogonal
#endif

/********************************** USAGE ************************************/
/*
//1. Header
#include <marching_cubes_33.h>
#include "MC33_util_grd.c"

//2. Read a grid file.
	_GRD* G = read_dat_file(stagbeetle832x832x494.dat);
	// The dataset stagbeetle832x832x494.dat can be downloaded at
	// https://www.cg.tuwien.ac.at/research/publications/2005/dataset-stagbeetle/
	// see MC33_util_grd.c for other grid loading functions

//3. Create a MC33 structure.
	MC33 *M = create_MC33(G);

//4. calculate an isosurface with isovalue v (a float).
	surface *S = calculate_isosurface(M, v);

//5. When finished calculating isosurfaces, free the memory occupied by M.
	free_MC33(M);
*/

/******************************* CUSTOMIZING *********************************/
//Do not change the following lines after compiling the library:

//#define integer_GRD // for dataset with integer type
//#define size_type_GRD 4 // 1, 2, 4 or 8 (8 for double, if not defined integer_GRD)

//#define GRD_orthogonal // If defined, the library only works with orthogonal grids.
/*****************************************************************************/
#define MC33C_VERSION_MAJOR 5
#define MC33C_VERSION_MINOR 3

#if defined(integer_GRD)
#if size_type_GRD == 4
/*
GRD_data_type is the variable type of the grid data, by default it is float.
*/
typedef unsigned int GRD_data_type;
#elif size_type_GRD == 2
typedef unsigned short int GRD_data_type;
#elif size_type_GRD == 1
typedef unsigned char GRD_data_type;
#else
#error "Incorrect size of the data type. size_type_GRD permitted values: 1, 2 or 4."
typedef float GRD_data_type;
#endif
#elif size_type_GRD == 8
typedef double GRD_data_type;
#else
typedef float GRD_data_type;
#undef size_type_GRD
#define size_type_GRD 4
#endif

#if !defined(marching_cubes_33_c) && defined(__cplusplus)
	extern "C" {
#endif

/*
The structure _GRD contains a function F[][][] evaluated at a grid of regularly
spaced points. N[] is the number of intervals in each dimension, so the number
of points is N[] + 1, and the grid contain (N[2] + 1)*(N[1] + 1)*(N[0] + 1)
points. L[] is the grid size in each dimension. r0[] are the coordinates of
the first grid point. d[] is the distance between adjacent points in each
dimension (can be different for each dimension), nonortho has a value of 1 when
the grid is inclined else 0. _A and A_ are the matrices that transform from
inclined to orthogonal coordinates and vice versa, respectively. If the grid is
periodic (is infinitely repeated along each dimension) the flag periodic must
be 1 (for x) + 2 (for y) + 4 (for z), else 0.

In this program, if GRD_orthogonal is defined, then nonortho, _A and A_ can be
removed from this structure, and only work with orthogonal grids. periodic
and L[] are not used here and also can be removed.
*/

#ifndef LAGRANGE3D4GRD_H
typedef struct
{
	GRD_data_type ***F;
	unsigned int N[3];
	double r0[3], d[3];
	float L[3]; 
#ifndef GRD_orthogonal
	float Ang[3];//not necessary
	int nonortho; //necessary if GRD_orthogonal is not defined
	double _A[3][3], A_[3][3]; //necessary if GRD_orthogonal is not defined
#endif
	int periodic;//not necessary
	int internal_data;
	char title[160];//not necessary
} _GRD;
#endif

/*The struct surface contains the data of an isosurface.
The number of points is nV, and the number of triangles is nT. The coordinates
of the points are stored in the array V[][3].
N[][3] is the array that contains the normal vectors calculated at the points
of the surface. The array color[] contains the color of each point.
T[][3] is the array that contains the sets of three point (vertex) indices
that form each triangle.*/
typedef struct
{
/* A is a mask (n & A is equivalent to n % dim1).
n >> E is equivalent to n / dim1. */
	unsigned int (*T)[3];
	float (*V)[3];
	float (*N)[3];
	int *color;
	unsigned int nV, nT;
	union {
		long long unsigned int ul[2];
		int i[4];
		char c[16];
		float f[4];
		double df[2];
	} user; // user data
} surface;


typedef struct
{
// copy of some variables of the surface structure:
	unsigned int (*T)[3];
	float (*V)[3];
	float (*N)[3];
	int *color;
	unsigned int nV, nT;
	float iso;

/*memoryfault takes a non-zero value if the system memory is insufficient
to store the surface.*/
	int memoryfault;

	unsigned int capt, capv;

// copy of some variables of the _GRD structure:
	const GRD_data_type ***F;
	float O[3], D[3], ca, cb;
	//int cubicflag;
	unsigned int nx, ny, nz;
	unsigned int (*store)(void *, float *);
#ifndef GRD_orthogonal
	double _A[3][3], A_[3][3];
#endif

// temporary structures that store the indices of triangle vertices:
	unsigned int **Dx, **Dy, **Ux, **Uy, **Lz;
} MC33;

#ifndef compiling_libMC33
/*
DefaultColorMC = 0x00FF0080
										 B G R
Red 128, green 0, blue 255.
*/
int DefaultColorMC = 0xff5c5c5c;//gray RGB color as unsigned char [3]

#ifndef GRD_orthogonal
#ifndef multTA_bf_code
#define multTA_bf_code
//c = Ab, A is a 3x3 upper triangular matrix. If t != 0, A is transposed.
void _multTSA_bf(const double (*A)[3], float *b, float *c, int t)
{
	if(t)
	{
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
		c[1] = A[0][1]*b[0] + A[1][1]*b[1];
		c[0] = A[0][0]*b[0];
	}
	else
	{
		c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		c[1] = A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][2]*b[2];
	}
}
//Performs the multiplication of the matrix A and the vector b: c = Ab. If t != 0, A is transposed.
void _multA_bf(const double (*A)[3], float* b, float* c, int t)
{
	double u,v;
	if(t)
	{
		u = A[0][0]*b[0] + A[1][0]*b[1] + A[2][0]*b[2];
		v = A[0][1]*b[0] + A[1][1]*b[1] + A[2][1]*b[2];
		c[2] = A[0][2]*b[0] + A[1][2]*b[1] + A[2][2]*b[2];
	}
	else
	{
		u = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
		v = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
		c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
	}
	c[0] = u;
	c[1] = v;
}
void (*mult_Abf)(const double (*)[3], float *, float *, int) = _multA_bf;
#endif // multTA_bf_code
#endif // GRD_orthogonal
#endif // compiling_libMC33

/******************************************************************
Saves all the surface *S data (in binary format) to a "filename" file. The
return value is 0 if the call succeeds, else -1.
*/
int write_bin_s(const char *filename, surface *S);
/******************************************************************
Reads (from a "filename" file) the surface data stored in binary format.
The return value is a surface pointer if the call succeeds, else NULL.
*/
surface* read_bin_s(const char *filename);
/******************************************************************
Saves all the surface *S data (in plain text format) to a "filename" file.
The return value is 0 if the call succeeds, else -1.
*/
int write_txt_s(const char *filename, surface *S);

/******************************************************************
Creates a MC33 structure, 
*/
MC33 *create_MC33(_GRD *G);

/******************************************************************
Calculates the isosurface (iso is the isovalue) using a MC33 structure
pointed by M. The return value is a pointer to surface. The pointer will
be NULL if there is not enough memory. The isovalue is stored in the
member user.f[3] of surface struct.*/
surface* calculate_isosurface(MC33 *M, float iso);

/******************************************************************
Return the size of surface without calculate it. The function can
calculate the number of vertices and triangles. If only the size is
required:
unsigned long long size = size_of_isosurface(M, v, 0, 0);
*/
unsigned long long size_of_isosurface(MC33 *M, float iso, unsigned int *nV, unsigned int *nT);

/******************************************************************
Releases the allocated memory occupied by MC33 structure pointed by M.
*/
void free_MC33(MC33 *M);

/******************************************************************
Releases the allocated memory pointed to by S.
*/
void free_surface_memory(surface *S);

#if !defined(marching_cubes_33_c) && defined(__cplusplus)
	} // extern "C"
#endif

#endif //marching_cubes_33_h

