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
	February 2026
	This library is a modified version of the library described in:
	Vega, D., Abache, J., Coll, D., A Fast and Memory Saving Marching Cubes 33
	implementation with the correct interior test, Journal of Computer Graphics
	Techniques (JCGT), vol. 8, no. 3, 1–18, 2019
*/

#ifndef marching_cubes_33_h
#define marching_cubes_33_h

/********************************** USAGE ************************************/
/*
//1. Header
#include <marching_cubes_33.h>
#include "MC33_util_grd.c"

//2. Read a grid file.
	_GRD* G = read_dat_file(stagbeetle832x832x494.dat);
	// The dataset stagbeetle832x832x494.dat can be downloaded at
	// https://www.cg.tuwien.ac.at/research/publications/2005/dataset-stagbeetle/
	// see MC33_util_grd.c for other grid reading functions

//3. Create a MC33 structure.
	MC33 *M = create_MC33(G);

//4. calculate an isosurface with isovalue v (a float).
	surface *S = calculate_isosurface(M, v);

//5. When finished, free the memory occupied by S, M, and G.
	free_surface_memory(S);
	free_MC33(M);
	free_memory_grd(G);
*/

/******************************* CUSTOMIZING *********************************/
//Do not change the following lines after compiling the library:

//#define INTEGER_GRD // Uncomment this define for dataset with integer type
//#define GRD_TYPE_SIZE 4 // 1, 2, 4 or 8 (8 for double, if not defined INTEGER_GRD)
//#define GRD_ORTHOGONAL // If defined, the library only works with orthogonal grids.
//#define MC_NORMAL_NEG // the front and back surfaces are exchanged.
//#define DEFAULT_SURFACE_COLOR 0xFF80FF40// RGBA 0xAABBGGRR: red 64, green 255, blue 128
/*****************************************************************************/

#define MC33C_VERSION_MAJOR 5
#define MC33C_VERSION_MINOR 4

#if defined(INTEGER_GRD)
#if GRD_TYPE_SIZE == 4
/*
GRD_data_type is the variable type of the grid data, by default it is float.
*/
typedef unsigned int GRD_data_type;
#elif GRD_TYPE_SIZE == 2
typedef unsigned short int GRD_data_type;
#elif GRD_TYPE_SIZE == 1
typedef unsigned char GRD_data_type;
#else
#error "Incorrect size of the data type. GRD_TYPE_SIZE permitted values: 1, 2 or 4."
typedef float GRD_data_type;
#endif
#elif GRD_TYPE_SIZE == 8
typedef double GRD_data_type;
#else
typedef float GRD_data_type;
#undef GRD_TYPE_SIZE
#define GRD_TYPE_SIZE 4
#endif

#if !defined(marching_cubes_33_c) && defined(__cplusplus) && !defined(mc33_no_lib)
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

In this program, if GRD_ORTHOGONAL is defined, then nonortho, _A and A_ can be
removed from this structure, and only work with orthogonal grids. periodic
and L[] are not used here and also can be removed.
*/

typedef struct {
	GRD_data_type ***F;
	unsigned int N[3];
	double r0[3], d[3];
	float L[3]; 
#ifndef GRD_ORTHOGONAL
	float Ang[3];//not necessary
	int nonortho; //necessary if GRD_ORTHOGONAL is not defined
	double _A[3][3], A_[3][3]; //necessary if GRD_ORTHOGONAL is not defined
#endif
	int periodic;//not necessary
	int internal_data;
	char title[160];//not necessary
} _GRD;

/*The struct surface contains the data of an isosurface.
The number of points is nV, and the number of triangles is nT. The coordinates
of the points are stored in the array V[][3].
N[][3] is the array that contains the normal vectors calculated at the points
of the surface. The array color[] contains the color of each point.
T[][3] is the array that contains the sets of three point (vertex) indices
that form each triangle.*/
typedef struct {
/* A is a mask (n & A is equivalent to n % dim1).
n >> E is equivalent to n / dim1. */
	unsigned int (*T)[3];
	float (*V)[3];
	float (*N)[3];
	int *color;
	unsigned int nV, nT;
	float iso;
	union {
		void *p;
		long long ul;
		int i[2];
		short si[4];
		char c[8];
		float f[2];
		double df;
	} user; // user data
} surface;

typedef struct {
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
	unsigned int nx, ny, nz;
	unsigned int (*store)(void *, float *);
#ifndef GRD_ORTHOGONAL
	double _A[3][3], A_[3][3];
#endif

// temporary structures that store the indices of triangle vertices:
	unsigned int **Dx, **Dy, **Ux, **Uy, **Lz;
} MC33;

extern int DefaultColorMC;
/* Example:
DefaultColorMC = 0x00FF0080
										 B G R
Red 128, green 0, blue 255.
*/

#ifndef GRD_ORTHOGONAL
extern void (*mult_Abf)(const double (*)[3], float *, float *, int);
#endif /* GRD_ORTHOGONAL */

/******************************************************************
Saves all the surface *S data (in binary format) to a "filename" file. The
return value is 0 if the call succeeds, else -1.
*/
int write_bin_s(surface *S, const char *filename);

/******************************************************************
Reads (from a "filename" file) the surface data stored in binary format.
The return value is a surface pointer if the call succeeds, else NULL.
*/
surface* read_bin_s(const char *filename);

/******************************************************************
Saves all the surface *S data (in plain text format) to a "filename" file.
The return value is 0 if the call succeeds, else -1.
*/
int write_txt_s(surface *S, const char *filename);

/******************************************************************
Saves the surface *S data (without the color) to Wavefront .obj file.
The return value is 0 if the call succeeds, else -1.
*/
int write_obj_s(surface *S, const char *filename);

/******************************************************************
Saves the surface *S data to Polygon File Format .ply file.
https://paulbourke.net/dataformats/ply/
The return value is 0 if the call succeeds, else -1.*/
int write_ply_s(surface *S, const char *filename, const char* author, const char* object);

/******************************************************************
Creates an MC33 structure from a pointer to a _GRD structure. If successful,
returns a pointer to the MC33 structure. On error, returns null pointer.
*/
MC33 *create_MC33(_GRD *G);

/******************************************************************
Calculates the isosurface (iso is the isovalue) using a MC33 structure
pointed by M. The return value is a pointer to surface. The pointer will
be NULL if there is not enough memory. The isovalue is stored in the
member iso of surface struct.*/
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

/******************************************************************
Free the memory occupied by the _GRD structure pointed to by Z.
*/
void free_memory_grd(_GRD *Z);

/******************************************************************
Allocate memory for the grid data. Before calling this function, memory must
be allocated for the _GRD structure and values must be assigned to the number
of points in each dimension (array N of the _GRD structure). The function
returns -1 if the allocation fails. Example:
	_GRD *Z = (_GRD *)malloc(sizeof(_GRD));
	Z->N[0] = 200; Z->N[1] = 200; Z->N[2] = 100;
	if (alloc_F(Z))
		return; // memory error
// fill Z->F[][][] here.
*/
int alloc_F(_GRD* Z);

/******************************************************************
read_grd read a filename file (the file must be a output *.grd file from the
DMol program), it returns a pointer to struct _GRD that contains all the grid
data.
*/
_GRD* read_grd(const char *filename);

/* Internal binary format
*/
_GRD* read_grd_binary(const char* filename);

/******************************************************************
Reads a set of files that contain a slab of res*res scan data points, the data
points are read as unsigned short int (if order is different from 0, the bytes
of the unsigned short are exchanged). The filename must end with a number, and
the fuction read all files with end number greater or equal to filename.
(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)
*/
_GRD* read_scanfiles(const char *filename, unsigned int res, int order);

/******************************************************************
Reads a file that contains only the data points as integers of 8, 16 or 32 bits.
byte is the number of bytes of each data point (1 to 4). If the data is big endian,
byte must be negative. The vector N[3] contains the number of points in each
dimension. The size of file must be byte*N[0]*N[1]*N[2].
*/
_GRD* read_raw_file(const char *filename, unsigned int *N, int byte, int isfloat);

/******************************************************************
Reads a dat file:
http://www.cg.tuwien.ac.at/research/vis/datasets/
*/
_GRD* read_dat_file(const char *filename);

/******************************************************************
	set_data_pointer creates a _GRD struct from an external data array. data
	must be stored with the nested inner loop running from i = 0 to Nx - 1
	and the outer loop from k = 0 to Nz - 1. The data will not be erased by
	free_memory_grd function. The function returns a pointer to the created
	struct.
*/
_GRD* grid_from_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data);

/******************************************************************
	Build a grid by using a scalar function
	double fn(double x, double y, double z)
*/
_GRD* generate_grid_from_fn(
	double x_initial, double y_initial, double z_initial,
	double x_final, double y_final, double z_final,
	double x_step, double y_step, double z_step,
	double (*fn)(double x, double y, double z));

#if !defined(marching_cubes_33_c) && defined(__cplusplus) && !defined(mc33_no_lib)
	} // extern "C"
#endif

#endif //marching_cubes_33_h

