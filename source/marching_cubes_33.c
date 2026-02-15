/*
	File: marching_cubes_33.c
	Programmed by: David Vega - dvega@uc.edu.ve
	               Javier Abache - jabache@uc.edu.ve
	March 2012
	updated by David Vega
	December 2014,
	June 2015,
	November 2017.
	February 2018.
	June 2019
	January 2020
	February 2020
	July 2020
	August 2020
	February 2021.
	August 2021
	February 2022
	February 2026
*/

#ifndef marching_cubes_33_c
#define marching_cubes_33_c

#if defined(marching_cubes_33_h) && !defined(compiling_libMC33) && defined(__cplusplus)
#error ***Do not include the file marching_cubes_33.h***
#include "*" //to abort the compilation
#endif

// Silence type conversion warnings in Visual C
#ifdef _MSC_VER
#pragma warning( push )
// disable warning when a double value is assigned to float variable
#pragma warning( disable : 4244 )
// disable warning when a size_t value is assigned to int variable
#pragma warning( disable : 4267 )
#endif

#include "MC33_LookUpTable.h"
#include "../include/marching_cubes_33.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#if defined (__SSE__) || ((defined (_M_IX86) || defined (_M_X64)) && !defined(_CHPE_ONLY_))
// https://stackoverflow.com/questions/59644197/inverse-square-root-intrinsics
// faster than 1.0f/std::sqrt, but with little accuracy.
#include <immintrin.h>
inline float invSqrt(float f) {
	__m128 temp = _mm_set_ss(f);
	temp = _mm_rsqrt_ss(temp);
	return _mm_cvtss_f32(temp);
}
#else
#include <math.h>
inline float invSqrt(float f) {
	return 1.0/sqrt(f);
}
#endif

#ifndef DEFAULT_SURFACE_COLOR
#define DEFAULT_SURFACE_COLOR 0xff5c5c5c;
#endif
int DefaultColorMC = DEFAULT_SURFACE_COLOR;

/****************** Surface managing functions ****************/

void free_surface_memory(surface *S) {
	if (S) {
		free(S->T);
		free(S->V);
		free(S->N);
		free(S->color);
		free(S);
	}
}

int write_bin_s(surface *S, const char *filename) {
	int i;
	FILE *out = fopen(filename,"wb");
	if (!out)
		return -1;
	fputs(".sup",out);

	fwrite(S->user.f + 3,sizeof(float),1,out);
	fwrite(&S->nV,sizeof(int),1,out);
	fwrite(&S->nT,sizeof(int),1,out);

	fwrite(S->T,3*S->nT*sizeof(int),1,out);
	fwrite(S->V,3*S->nV*sizeof(float),1,out);
	fwrite(S->N,3*S->nV*sizeof(float),1,out);
	i = fwrite(S->color,S->nV*sizeof(int),1,out);
	fclose(out);
	return -(i != 1);
}

surface* read_bin_s(const char *filename) {
	surface *S;
	int i;

	FILE *in = fopen(filename,"rb");
	if (!in)
		return 0;
	fread(&i,sizeof(int),1,in);
	S = (surface*)malloc(sizeof(surface));
	if (i != 0x7075732e || !S) {
		fclose(in);
		free(S);
		return 0;
	}
	fread(S->user.f + 3,sizeof(float),1,in);
	fread(&S->nV,sizeof(int),1,in);
	fread(&S->nT,sizeof(int),1,in);

	S->T = (unsigned int (*)[3])malloc(S->nT*sizeof(void*));
	S->V = (float (*)[3])malloc(S->nV*sizeof(void*));
	S->N = (float (*)[3])malloc(S->nV*sizeof(void*));
	S->color = (int *)malloc(S->nV*sizeof(void*));
	if (!(S->T && S->V && S->N && S->color)) {
		free_surface_memory(S);
		fclose(in);
		return 0;
	}
	i = !fread(S->T,3*S->nT*sizeof(int),1,in);
	i += !fread(S->V,3*S->nV*sizeof(float),1,in);
	i += !fread(S->N,3*S->nV*sizeof(float),1,in);
	i += !fread(S->color,S->nV*sizeof(int),1,in);
	fclose(in);
	if (i) {
		free_surface_memory(S);
		return 0;
	}
	return S;
}

int write_txt_s(surface *S, const char *filename) {
	FILE *out;
	unsigned int i, *t;
	float *r;

	out = fopen(filename,"w");
	if (!out)
		return -1;

	fprintf(out,"isovalue: %10.5E\n\nVERTICES:\n",S->iso);
	fprintf(out,"%d\n\n",S->nV);
	for (i = 0; i != S->nV; i++) {
		r = S->V[i];
		fprintf(out,"%9.6f %9.6f %9.6f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"\n\nTRIANGLES:\n");
	fprintf(out,"%d\n\n",S->nT);
	for (i = 0; i != S->nT; i++) {
		t = S->T[i];
		fprintf(out,"%8d %8d %8d\n",t[0],t[1],t[2]);
	}

	fprintf(out,"\n\nNORMALS:\n");
	for (i = 0; i != S->nV; i++) {
		r = S->N[i];
		fprintf(out,"%9.6f %9.6f %9.6f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"\n\nCOLORS:\n");
	for (i = 0; i != S->nV; i++)
		fprintf(out,"%d\n",S->color[i]);
	i = fprintf(out,"\nEND\n");
	fclose(out);
	return -(i < 5);
}

int write_obj_s(surface *S, const char *filename) {
	FILE *out;
	unsigned int i, *t;
	float *r;
	char s0[12], s1[12], s2[12];

	out = fopen(filename,"w");
	if (!out)
		return -1;

	fprintf(out,"# isovalue: %10.5E\n# VERTICES %d:\n", S->iso, S->nV);
	for (i = 0; i != S->nV; i++) {
		r = S->V[i];
		fprintf(out,"v %f %f %f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"# NORMALS:\n");
	for (i = 0; i != S->nV; i++) {
		r = S->N[i];
		fprintf(out,"vn %f %f %f\n",r[0],r[1],r[2]);
	}

	fprintf(out,"# TRIANGLES %d:\n",S->nT);
	for (i = 0; i != S->nT; i++) {
	    t = S->T[i];
		sprintf(s0,"%d", t[0] + 1);
		sprintf(s1,"%d", t[1] + 1);
		sprintf(s2,"%d", t[2] + 1);
		fprintf(out,"f %s//%s %s//%s %s//%s\n", s0, s0, s1, s1, s2, s2);
	}

	i = fprintf(out,"# END");
	fclose(out);
	return -(i < 5);
}

int write_ply_s(surface *S, const char *filename, const char* author, const char* object) {
	FILE *out;
	unsigned int i, *t;
	float *r;
	unsigned char* c;
	char empty = 0;

	out = fopen(filename,"w");
	if (!out)
		return -1;

	if (!author)
		author = &empty;
	if (!object)
		object = &empty;

	fprintf(out,"ply\nformat ascii 1.0\ncomment author: %s\ncomment object: %s\n",author,object);
	fprintf(out,"element vertex %d\nproperty float x\nproperty float y\nproperty float z",S->nV);
	fprintf(out,"\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face");
	fprintf(out," %d\nproperty list uchar int vertex_index\nend_header",S->nT);

	for (i = 0; i < S->nV; i++) {
		r = S->V[i];
		c = (unsigned char*)(S->color + i);
		fprintf(out,"\n%f %f %f %d %d %d",r[0],r[1],r[2],c[0],c[1],c[2]);
	}

	for (i = 0; i < S->nT; i++) {
		t = S->T[i];
		fprintf(out,"\n3 %d %d %d",t[0],t[1],t[2]);
	}

	i = fprintf(out,"\n");
	fclose(out);
	return -(!i);
}

/***************** Marching cubes 33 functions ****************/

/******************************************************************
Vertices:           Faces:
    3 __________2        ___________
   /|          /|      /|          /|
  / |         / |     / |   2     / |
7/__________6/  |    /  |     4  /  |
|   |       |   |   |¯¯¯¯¯¯¯¯¯¯¯| 1 |     z
|   0_______|___1   | 3 |_______|___|     |
|  /        |  /    |  /  5     |  /      |____y
| /         | /     | /     0   | /      /
4/__________5/      |/__________|/      x


This function return a vector with all six test face results (face[6]). Each
result value is 1 if the positive face vertices are joined, -1 if the negative
vertices are joined, and 0 (unchanged) if the test must no be applied. The
return value of this function is the the sum of all six results.*/
int MC33_faceTests(int *face, int ind, const float *v) {
	if (ind&0x80)//vertex 0
	{
		face[0] = ((ind&0xCC) == 0x84? (v[0]*v[5] < v[1]*v[4]? -1: 1): 0);//0x84 = 10000100, vertices 0 and 5
		face[3] = ((ind&0x99) == 0x81? (v[0]*v[7] < v[3]*v[4]? -1: 1): 0);//0x81 = 10000001, vertices 0 and 7
		face[4] = ((ind&0xF0) == 0xA0? (v[0]*v[2] < v[1]*v[3]? -1: 1): 0);//0xA0 = 10100000, vertices 0 and 2
	}
	else
	{
		face[0] = ((ind&0xCC) == 0x48? (v[0]*v[5] < v[1]*v[4]? 1: -1): 0);//0x48 = 01001000, vertices 1 and 4
		face[3] = ((ind&0x99) == 0x18? (v[0]*v[7] < v[3]*v[4]? 1: -1): 0);//0x18 = 00011000, vertices 3 and 4
		face[4] = ((ind&0xF0) == 0x50? (v[0]*v[2] < v[1]*v[3]? 1: -1): 0);//0x50 = 01010000, vertices 1 and 3
	}
	if (ind&0x02)//vertex 6
	{
		face[1] = ((ind&0x66) == 0x42? (v[1]*v[6] < v[2]*v[5]? -1: 1): 0);//0x42 = 01000010, vertices 1 and 6
		face[2] = ((ind&0x33) == 0x12? (v[3]*v[6] < v[2]*v[7]? -1: 1): 0);//0x12 = 00010010, vertices 3 and 6
		face[5] = ((ind&0x0F) == 0x0A? (v[4]*v[6] < v[5]*v[7]? -1: 1): 0);//0x0A = 00001010, vertices 4 and 6
	}
	else
	{
		face[1] = ((ind&0x66) == 0x24? (v[1]*v[6] < v[2]*v[5]? 1: -1): 0);//0x24 = 00100100, vertices 2 and 5
		face[2] = ((ind&0x33) == 0x21? (v[3]*v[6] < v[2]*v[7]? 1: -1): 0);//0x21 = 00100001, vertices 2 and 7
		face[5] = ((ind&0x0F) == 0x05? (v[4]*v[6] < v[5]*v[7]? 1: -1): 0);//0x05 = 00000101, vertices 5 and 7
	}
	return face[0] + face[1] + face[2] + face[3] + face[4] + face[5];
}

/* Faster function for the face test, the test is applied to only one face
(int face). This function is only used for the cases 3 and 6 of MC33*/
int MC33_faceTest1(int face, const float *v) {
	switch (face) {
	case 0:
		return (v[0]*v[5] < v[1]*v[4]? 0x48: 0x84);
	case 1:
		return (v[1]*v[6] < v[2]*v[5]? 0x24: 0x42);
	case 2:
		return (v[3]*v[6] < v[2]*v[7]? 0x21: 0x12);
	case 3:
		return (v[0]*v[7] < v[3]*v[4]? 0x18: 0x81);
	case 4:
		return (v[0]*v[2] < v[1]*v[3]? 0x50: 0xA0);
	default:
		return (v[4]*v[6] < v[5]*v[7]? 0x05: 0x0A);
	}
}

// an ugly signbit for float type with
// warning: dereferencing type-punned pointer will break strict-aliasing rules
// Silence dereferencing type-punned pointer warning in GCC
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

inline unsigned int signbf(float x) {
	return ((*(unsigned int*)(void*)(&x))&0x80000000);
}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/******************************************************************
Interior test function. If the test is positive, the function returns a value
different from 0. The integer i must be 0 to test if the vertices 0 and 6 are
joined. 1 for vertices 1 and 7, 2 for vertices 2 and 4, and 3 for 3 and 5.
For case 13, the integer flag13 must be 1, and the function returns 2 if one
of the vertices 0, 1, 2 or 3 is joined to the center point of the cube (case
13.5.2), returns 1 if one of the vertices 4, 5, 6 or 7 is joined to the
center point of the cube (case 13.5.2 too), and it returns 0 if the vertices
are no joined (case 13.5.1)*/
int MC33_interiorTest(int i, int flag13, const float *v) {
	//Signs of cube vertices were changed to use signbit function in calc_isosurface
	//A0 = -v[0], B0 = -v[1], C0 = -v[2], D0 = -v[3]
	//A1 = -v[4], B1 = -v[5], C1 = -v[6], D1 = -v[7]
	//But the function still works
	float At = v[4] - v[0], Bt = v[5] - v[1], Ct = v[6] - v[2], Dt = v[7] - v[3];
	float t = At*Ct - Bt*Dt; // the "a" value.
	if (signbf(t)) {
		if (i&0x01)
			return 0;
	} else {
		if (!(i&0x01) || t == 0)
			return 0;
	}
	t = 0.5f*(v[3]*Bt - v[2]*At + v[1]*Dt - v[0]*Ct)/t; // t = -b/2a
	if (t > 0 && t < 1) {
		At = v[0] + At*t;
		Bt = v[1] + Bt*t;
		Ct = v[2] + Ct*t;
		Dt = v[3] + Dt*t;
		Ct *= At;
		Dt *= Bt;
		if (i&0x01) {
			if (Ct < Dt && signbf(Dt) == 0)
				return (signbf(Bt) == signbf(v[i])) + flag13;
		} else {
			if (Ct > Dt && signbf(Ct) == 0)
				return (signbf(At) == signbf(v[i])) + flag13;
		}
	}
	return 0;
}

unsigned int MC33_fail_mem_VN(MC33 *M) {
	M->nV = 0;
	M->memoryfault = 1;
	return 0;
}

void MC33_fail_mem_T(MC33 *M) {
	M->nT = 0;
	M->capt >>= 1;
	M->memoryfault = 1;
	return;
}

/******************************************************************
Assign memory for the vertex r[3], normal n[3]. The return value is the new
vertex label.
*/
unsigned int MC33_spn0(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	for (int i = 0; i != 3; i++)
		p[i] = *(r++);
	// now r points to normal coordinates
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
unsigned int MC33_spnA(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	for (int i = 0; i != 3; i++)
		p[i] = *(r++)*M->D[i] + M->O[i];
	// now r points to normal coordinates
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
unsigned int MC33_spnB(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	for (int i = 0; i != 3; i++)
		p[i] = *(r++)*M->D[i] + M->O[i];
	// now r points to normal coordinates
	r[0] *= M->ca; // normal[0]
	r[1] *= M->cb; // normal[1]
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
#ifndef GRD_ORTHOGONAL
unsigned int MC33_spnC(void *mc33, float *r) {
	MC33 *M = (MC33 *)mc33;
	unsigned int nv = M->nV++;
	float t, *p;
	if (nv == M->capv) {
		void *pt;
		pt = realloc(M->N,nv*6*sizeof(float)); // the memory space is duplicated
		if (pt) {
			M->N = (float(*)[3])pt;
			pt = realloc(M->V,nv*6*sizeof(float));
			if (pt) {
				M->V = (float(*)[3])pt;
				M->capv <<= 1;
			}
			else
				return MC33_fail_mem_VN(M);
		} else
			return MC33_fail_mem_VN(M);
	}
	p = M->V[nv];
	mult_Abf(M->_A,r,r,0);
	for (int i = 0; i != 3; i++)
		p[i] = *(r++) + M->O[i];
	// now r points to normal coordinates
	mult_Abf(M->A_,r,r,1);
#ifndef MC_NORMAL_NEG
	t = invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#else //MC_NORMAL_NEG reverse the direction of the normal
	t = -invSqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
#endif
	p = M->N[nv];
	*p = t * *r; *(++p) = t * *(++r); *(++p) = t * *(++r);
	return nv;
}
#endif

/******************************************************************
Auxiliary function that calculates the normal if a vertex
of the cube lies on the isosurface.
*/
unsigned int MC33_surfint(MC33 *M, unsigned int x, unsigned int y, unsigned int z, float *r) {
	r[0] = x; r[1] = y; r[2] = z;
	if (x == 0)
		r[3] = M->F[z][y][0] - M->F[z][y][1];
	else if (x == M->nx)
		r[3] = M->F[z][y][x - 1] - M->F[z][y][x];
	else
		r[3] = 0.5f*(M->F[z][y][x - 1] - M->F[z][y][x + 1]);
	if (y == 0)
		r[4] = M->F[z][0][x] - M->F[z][1][x];
	else if (y == M->ny)
		r[4] = M->F[z][y - 1][x] - M->F[z][y][x];
	else
		r[4] = 0.5f*(M->F[z][y - 1][x] - M->F[z][y + 1][x]);
	if (z == 0)
		r[5] = M->F[0][y][x] - M->F[1][y][x];
	else if (z == M->nz)
		r[5] = M->F[z - 1][y][x] - M->F[z][y][x];
	else
		r[5] = 0.5f*(M->F[z - 1][y][x] - M->F[z + 1][y][x]);
	return M->store(M, r);
}

/******************************************************************
This function find the MC33 case (using the index i, and the face and interior
tests). The correct triangle pattern is obtained from the arrays contained in
the MC33_LookUpTable.h file. The necessary vertices (intersection points) are
also calculated here (using trilinear interpolation).
       _____2_____
     /|          /|
   11 |<-3     10 |
   /____6_____ /  1     z
  |   |       |   |     |
  |   |_____0_|___|     |____y
  7  /        5  /     /
  | 8         | 9     x
  |/____4_____|/

The temporary matrices: M->Lz, M->Dx, M->Dy, M->Ux and M->Uy are filled
and used here.*/
#define FF 0xFFFFFFFF
void MC33_findCase(MC33 *M, unsigned int x, unsigned int y, unsigned int z, unsigned int i, const float *v) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	unsigned int ti[3];//for vertex indices of a triangle
	union { // memory saving
		int f[6];//for the face tests
		float r[6];//for intercept and normal coordinates
	} u;
	const unsigned short int *pcase = MC33_all_tables;
	unsigned int c, m, k, n;
	float t;
	if (i&0x80) {
		c = pcase[i^0xFF];
		m = (c&0x800) == 0;
		n = !m;
	} else {
		c = pcase[i];
		n = (c&0x800) == 0;
		m = !n;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // case 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&MC33_faceTest1(k>>2, v)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (MC33_interiorTest(k, 0, v)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&MC33_faceTest1(k%6, v))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (MC33_interiorTest(k/6, 0, v)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (u.f[4] + u.f[5] < 0? (u.f[0] + u.f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (u.f[4] + u.f[5] < 0? 983: (u.f[0] + u.f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (MC33_interiorTest(k>>1, 0, v)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2:
				if (k == 2? MC33_interiorTest(0, 0, v): MC33_interiorTest(0, 0, v)||MC33_interiorTest(k? 1: 3, 0, v))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (u.f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k == 2? MC33_interiorTest(1, 0, v): MC33_interiorTest(2, 0, v)||MC33_interiorTest(k? 3: 1, 0, v))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2: //12.1
				pcase += (MC33_interiorTest((0xDA010C>>(k<<1))&3, 0, v)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (u.f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (MC33_interiorTest((0xA7B7E5>>(k<<1))&3, 0, v)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(MC33_faceTests(u.f, 165, v))) {
			case 0:
				k = ((u.f[1] < 0)<<1)|(u.f[5] < 0);
				if (u.f[0]*u.f[1] == u.f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = MC33_interiorTest(k, 1, v); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 + (c? 10*k - 40*c: 6*k);
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((u.f[0] < 0? u.f[2] > 0: 12 + (u.f[2] < 0)) + (u.f[1] < 0? u.f[3] < 0: 6 + (u.f[3] > 0)));
				if (u.f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*u.f[0] + 4*u.f[1] + 3*u.f[2] + 2*u.f[3] + u.f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*u.f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) { // the vertices r[3] and normals (r + 3)[3] are calculated here
				case 0:
					if (z || x)
						ti[--k] = p[0] = M->Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = M->Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = M->Dx[y][0];
							else if (y? signbf(M->iso - M->F[0][y - 1][0]): 0)
								p[0] = M->Dy[y - 1][0];
							else
								p[0] = MC33_surfint(M,0,y,0,u.r);
						} else if (v[1] == 0) {
							if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = (p[1] != FF? p[1]: MC33_surfint(M,0,y + 1,0,u.r));
						} else {
							t = v[0]/(v[0] - v[1]);
							u.r[0] = u.r[2] = 0;
							u.r[1] = y + t;
							u.r[3] = (v[4] - v[0])*(1 - t) + (v[5] - v[1])*t;
							u.r[4] = v[1] - v[0];
							u.r[5] = (v[3] - v[0])*(1 - t) + (v[2] - v[1])*t;
							p[0] = M->store(M, u.r);
						}
						M->Dy[y][0] = ti[--k] = p[0];
					}
					break;
				case 1:
					if (x)
						ti[--k] = p[1] = M->Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = M->Dy[y][0];
							//else if (z && signbf(v[5]))
							//	p[1] = M->Dx[y + 1][0];
							else if (z && y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][0]): 0)
								p[1] = M->Dy[y + 1][0];
							else if (z? signbf(M->iso - M->F[z - 1][y + 1][0]): 0) {
								ti[--k] = p[1] = M->Lz[y + 1][0]; // value of previous slice
								break;
							} else
								p[1] = MC33_surfint(M, 0, y + 1, z, u.r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = (p[2] != FF? p[2]: MC33_surfint(M, 0, y + 1, z + 1, u.r));
						} else {
							t = v[1]/(v[1] - v[2]);
							u.r[0] = 0; u.r[1] = y + 1;
							u.r[2] = z + t;
							u.r[3] = (v[5] - v[1])*(1 - t) + (v[6] - v[2])*t;
							u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[z][y][0] - M->F[z][y + 2][0])*(1 - t)
										+ (M->F[z + 1][y][0] - M->F[z + 1][y + 2][0])*t):
										(v[1] - v[0])*(1 - t) + (v[2] - v[3])*t);
							u.r[5] = v[2] - v[1];
							p[1] = M->store(M, u.r);
						}
						M->Lz[y + 1][0] = ti[--k] = p[1];
					}
					break;
				case 2:
					if (x)
						ti[--k] = p[2] = M->Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = M->Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = M->Ux[y][0];
							else if (y? signbf(M->iso - M->F[z + 1][y - 1][0]): 0)
								p[2] = M->Uy[y - 1][0];
							else
								p[2] = MC33_surfint(M,0,y,z + 1,u.r);
						} else if (v[2] == 0) {
							if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = (p[1] != FF? p[1]: MC33_surfint(M,0,y + 1,z + 1,u.r));
						} else {
							t = v[3]/(v[3] - v[2]);
							u.r[0] = 0; u.r[2] = z + 1;
							u.r[1] = y + t;
							u.r[3] = (v[7] - v[3])*(1 - t) + (v[6] - v[2])*t;
							u.r[4] = v[2] - v[3];
							u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][y][0] - M->F[z + 2][y][0])*(1 - t)
										+ (M->F[z][y + 1][0] - M->F[z + 2][y + 1][0])*t):
										(v[3] - v[0])*(1 - t) + (v[2] - v[1])*t);
							p[2] = M->store(M, u.r);
						}
						M->Uy[y][0] = ti[--k] = p[2];
					}
					break;
				case 3:
					if (y || x)
						ti[--k] = p[3] = M->Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = M->Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = M->Dx[0][0];
							else if (z? signbf(M->iso - M->F[z - 1][0][0]): 0) {
								ti[--k] = p[3] = M->Lz[0][0]; // value of previous slice
								break;
							} else
								p[3] = MC33_surfint(M,0,0,z,u.r);
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else
								p[3] = (p[11] != FF? p[11]: MC33_surfint(M,0,0,z + 1,u.r));
						} else {
							t = v[0]/(v[0] - v[3]);
							u.r[0] = u.r[1] = 0;
							u.r[2] = z + t;
							u.r[3] = (v[4] - v[0])*(1 - t) + (v[7] - v[3])*t;
							u.r[4] = (v[1] - v[0])*(1 - t) + (v[2] - v[3])*t;
							u.r[5] = v[3] - v[0];
							p[3] = M->store(M,u.r);
						}
						M->Lz[0][0] = ti[--k] = p[3];
					}
					break;
				case 4:
					if (z)
						ti[--k] = p[4] = M->Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[0]))
								p[4] = M->Dx[y][x];
							else if (y && signbf(v[7]))
								p[4] = M->Lz[y][x + 1];
							else if (y? signbf(M->iso - M->F[0][y - 1][x + 1]): 0)
								p[4] = M->Dy[y - 1][x + 1];
							else if (y && x + 1 < M->nx? signbf(M->iso - M->F[0][y][x + 2]): 0)
								p[4] = M->Dx[y][x + 1];
							else
								p[4] = MC33_surfint(M,x + 1,y,0,u.r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[4] = p[5];
							else
								p[4] = (p[9] != FF? p[9]: MC33_surfint(M,x + 1,y + 1,0,u.r));
						} else {
							t = v[4]/(v[4] - v[5]);
							u.r[0] = x + 1; u.r[2] = 0;
							u.r[1] = y + t;
							u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[0][y][x] - M->F[0][y][x + 2])*(1 - t)
										+ (M->F[0][y + 1][x] - M->F[0][y + 1][x + 2])*t):
										(v[4] - v[0])*(1 - t) + (v[5] - v[1])*t);
							u.r[4] = v[5] - v[4];
							u.r[5] = (v[7] - v[4])*(1 - t) + (v[6] - v[5])*t;
							p[4] = M->store(M,u.r);
						}
						M->Dy[y][x + 1] = ti[--k] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (z) {
							if (signbf(v[4]))
								p[5] = p[4] = M->Dy[y][x + 1];
							else if (signbf(v[1]))
								p[5] = p[9] = M->Dx[y + 1][x];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z][y + 1][x + 2]): 0)
								p[5] = M->Dx[y + 1][x + 1];
							else if (y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][x + 1]): 0)
								p[5] = M->Dy[y + 1][x + 1];
							else if (signbf(M->iso - M->F[z - 1][y + 1][x + 1])) {
								ti[--k] = p[5] = M->Lz[y + 1][x + 1]; // value of previous slice
								break;
							} else
								p[5] = MC33_surfint(M,x + 1,y + 1,z,u.r);
						} else
							p[5] = MC33_surfint(M,x + 1,y + 1,0,u.r);
					} else if (v[6] == 0)
						p[5] = MC33_surfint(M,x + 1,y + 1,z + 1,u.r);
					else {
						t = v[5]/(v[5] - v[6]);
						u.r[0] = x + 1; u.r[1] = y + 1;
						u.r[2] = z + t;
						u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[z][y + 1][x] - M->F[z][y + 1][x + 2])*(1 - t)
									+ (M->F[z + 1][y + 1][x] - M->F[z + 1][y + 1][x + 2])*t):
									(v[5] - v[1])*(1 - t) + (v[6] - v[2])*t);
						u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[z][y][x + 1] - M->F[z][y + 2][x + 1])*(1 - t)
									+ (M->F[z + 1][y][x + 1] - M->F[z + 1][y + 2][x + 1])*t):
									(v[5] - v[4])*(1 - t) + (v[6] - v[7])*t);
						u.r[5] = v[6] - v[5];
						p[5] = M->store(M,u.r);
					}
					M->Lz[y + 1][x + 1] = ti[--k] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = M->Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = M->Lz[y][x + 1];
							else if (signbf(M->iso - M->F[z + 1][y - 1][x + 1]))
								p[6] = M->Uy[y - 1][x + 1];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z + 1][y][x + 2]): 0)
								p[6] = M->Ux[y][x + 1];
							else
								p[6] = MC33_surfint(M,x + 1,y,z + 1,u.r);
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = MC33_surfint(M,x + 1,0,z + 1,u.r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[6] = p[5];
						else
							p[6] = (p[10] == FF? MC33_surfint(M,x + 1,y + 1,z + 1,u.r): p[10]);
					} else {
						t = v[7]/(v[7] - v[6]);
						u.r[0] = x + 1;
						u.r[1] = y + t; u.r[2] = z + 1;
						u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[z + 1][y][x] - M->F[z + 1][y][x + 2])*(1 - t)
									+ (M->F[z + 1][y + 1][x] - M->F[z + 1][y + 1][x + 2])*t):
									(v[7] - v[3])*(1 - t) + (v[6] - v[2])*t);
						u.r[4] = v[6] - v[7];
						u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][y][x + 1] - M->F[z + 2][y][x + 1])*(1 - t)
										+ (M->F[z][y + 1][x + 1] - M->F[z + 2][y + 1][x + 1])*t):
						(v[7] - v[4])*(1 - t) + (v[6] - v[5])*t);
						p[6] = M->store(M,u.r);
					}
					M->Uy[y][x + 1] = ti[--k] = p[6];
					break;
				case 7:
					if (y)
						ti[--k] = p[7] = M->Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = M->Dx[0][x];
							else if (z && signbf(v[5]))
								p[7] = M->Dy[0][x + 1];
							else if (z && x + 1 < M->nx? signbf(M->iso - M->F[z][0][x + 2]): 0)
								p[7] = M->Dx[0][x + 1];
							else if (z? signbf(M->iso - M->F[z - 1][0][x + 1]): 0) {
								ti[--k] = p[7] = M->Lz[0][x + 1]; // value of previous slice
								break;
							} else
								p[7] = MC33_surfint(M,x + 1,0,z,u.r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else
								p[7] = (p[11] != FF? p[11]: MC33_surfint(M,x + 1,0,z + 1,u.r));
						} else {
							t = v[4]/(v[4] - v[7]);
							u.r[0] = x + 1; u.r[1] = 0;
							u.r[2] = z + t;
							u.r[3] = (x + 1 < M->nx? 0.5f*((M->F[z][0][x] - M->F[z][0][x + 2])*(1 - t)
										+ (M->F[z + 1][0][x] - M->F[z + 1][0][x + 2])*t):
										(v[4] - v[0])*(1 - t) + (v[7] - v[3])*t);
							u.r[4] = (v[5] - v[4])*(1 - t) + (v[6] - v[7])*t;
							u.r[5] = v[7] - v[4];
							p[7] = M->store(M,u.r);
						}
						M->Lz[0][x + 1] = ti[--k] = p[7];
					}
					break;
				case 8:
					if (z || y)
						ti[--k] = p[8] = M->Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = M->Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = M->Dy[0][x];
							else if (x? signbf(M->iso - M->F[0][0][x - 1]): 0)
								p[8] = M->Dx[0][x - 1];
							else
								p[8] = MC33_surfint(M,x,0,0,u.r);
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else
								p[8] = (p[7] != FF? p[7]: MC33_surfint(M,x + 1,0,0,u.r));
						} else {
							t = v[0]/(v[0] - v[4]);
							u.r[1] = u.r[2] = 0;
							u.r[0] = x + t;
							u.r[3] = v[4] - v[0];
							u.r[4] = (v[1] - v[0])*(1 - t) + (v[5] - v[4])*t;
							u.r[5] = (v[3] - v[0])*(1 - t) + (v[7] - v[4])*t;
							p[8] = M->store(M,u.r);
						}
						M->Dx[0][x] = ti[--k] = p[8];
					}
					break;
				case 9:
					if (z)
						ti[--k] = p[9] = M->Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							//else if (p[1] != FF)
							//	p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = M->Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = M->Lz[y + 1][x];
							else if (x? signbf(M->iso - M->F[0][y + 1][x - 1]): 0)
								p[9] = M->Dx[y + 1][x - 1];
							else
								p[9] = MC33_surfint(M,x,y + 1,0,u.r);
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else
								p[9] = (p[4] != FF? p[4]: MC33_surfint(M,x + 1,y + 1,0,u.r));
						} else {
							t = v[1]/(v[1] - v[5]);
							u.r[1] = y + 1; u.r[2] = 0;
							u.r[0] = x + t;
							u.r[3] = v[5] - v[1];
							u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[0][y][x] - M->F[0][y + 2][x])*(1 - t)
										+ (M->F[0][y][x + 1] - M->F[0][y + 2][x + 1])*t):
										(v[1] - v[0])*(1 - t) + (v[5] - v[4])*t);
							u.r[5] = (v[2] - v[1])*(1 - t) + (v[6] - v[5])*t;
							p[9] = M->store(M,u.r);
						}
						M->Dx[y + 1][x] = ti[--k] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = M->Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = M->Uy[y][x];
							else if (signbf(M->iso - M->F[z + 1][y + 1][x - 1]))
								p[10] = M->Ux[y + 1][x - 1];
							else
								p[10] = MC33_surfint(M,x,y + 1,z + 1,u.r);
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = MC33_surfint(M,0,y + 1,z + 1,u.r);
					} else if (v[6] == 0) {
						if (p[5] != FF)
							p[10] = p[5];
						else
							p[10] = (p[6] != FF? p[6]: MC33_surfint(M,x + 1,y + 1,z + 1,u.r));
					} else {
						t = v[2]/(v[2] - v[6]);
						u.r[0] = x + t;
						u.r[1] = y + 1; u.r[2] = z + 1;
						u.r[3] = v[6] - v[2];
						u.r[4] = (y + 1 < M->ny? 0.5f*((M->F[z + 1][y][x] - M->F[z + 1][y + 2][x])*(1 - t)
									+ (M->F[z + 1][y][x + 1] - M->F[z + 1][y + 2][x + 1])*t):
									(v[2] - v[3])*(1 - t) + (v[6] - v[7])*t);
						u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][y + 1][x] - M->F[z + 2][y + 1][x])*(1 - t)
									+ (M->F[z][y + 1][x + 1] - M->F[z + 2][y + 1][x + 1])*t):
									(v[2] - v[1])*(1 - t) + (v[6] - v[5])*t);
						p[10] = M->store(M,u.r);
					}
					M->Ux[y + 1][x] = ti[--k] = p[10];
					break;
				case 11:
					if (y)
						ti[--k] = p[11] = M->Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = M->Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = M->Uy[0][x];
							else if (x? signbf(M->iso - M->F[z + 1][0][x - 1]): 0)
								p[11] = M->Ux[0][x - 1];
							else
								p[11] = MC33_surfint(M,x,0,z + 1,u.r);
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else
								p[11] = (p[7] != FF? p[7]: MC33_surfint(M,x + 1,0,z + 1,u.r));
						} else {
							t = v[3]/(v[3] - v[7]);
							u.r[1] = 0; u.r[2] = z + 1;
							u.r[0] = x + t;
							u.r[3] = v[7] - v[3];
							u.r[4] = (v[2] - v[3])*(1 - t) + (v[6] - v[7])*t;
							u.r[5] = (z + 1 < M->nz? 0.5f*((M->F[z][0][x] - M->F[z + 2][0][x])*(1 - t)
										+ (M->F[z][0][x + 1] - M->F[z + 2][0][x + 1])*t):
										(v[3] - v[0])*(1 - t) + (v[7] - v[4])*t);
							p[11] = M->store(M,u.r);
						}
						M->Ux[0][x] = ti[--k] = p[11];
					}
				break;
				default:
					u.r[0] = x + 0.5f; u.r[1] = y + 0.5f; u.r[2] = z + 0.5f;
					u.r[3] = v[4] + v[5] + v[6] + v[7] - v[0] - v[1] - v[2] - v[3];
					u.r[4] = v[1] + v[2] + v[5] + v[6] - v[0] - v[3] - v[4] - v[7];
					u.r[5] = v[2] + v[3] + v[6] + v[7] - v[0] - v[1] - v[4] - v[5];
					ti[--k] = p[12] = M->store(M,u.r);
				}
			} else
				ti[--k] = p[c];//now ti contains the vertex indices of the triangle
		}
		if (ti[0] != ti[1] && ti[0] != ti[2] && ti[1] != ti[2]) { //to avoid zero area triangles
			if (M->nT == M->capt) {
				unsigned int (*pt)[3] = M->T;
				M->capt <<= 1;
				M->T = (unsigned int (*)[3])realloc(pt,M->capt*3*sizeof(int));
				if (!M->T) {
					M->T = pt;
					MC33_fail_mem_T(M);
				}
			}
			unsigned int *vp = M->T[M->nT++];
#ifndef MC_NORMAL_NEG
			*vp = ti[n]; *(++vp) = ti[m]; *(++vp) = ti[2];
#else
			*vp = ti[m]; *(++vp) = ti[n]; *(++vp) = ti[2];
#endif
		}
	}
}

void MC33_Case_count(MC33 *M, unsigned int x, unsigned int y, unsigned int z, unsigned int i, const float *v) {
	unsigned int p[13] = {FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF,FF};
	union {
		int f[6];
		unsigned int ti[3];
	} u;
	const unsigned short int *pcase = MC33_all_tables;
	unsigned int c, m, k;
	if (i&0x80) {
		c = pcase[i^0xFF];
		m = (c&0x800) == 0;
	} else {
		c = pcase[i];
		m = (c&0x800) != 0;
	}
	k = c&0x7FF;
	switch (c>>12) { //find the MC33 case
		case 0: // cases 1, 2, 5, 8, 9, 11 and 14
			pcase += k;
			break;
		case 1: // case 3
			pcase += ((m? i: i^0xFF)&MC33_faceTest1(k>>2, v)? 183 + (k<<1): 159 + k);
			break;
		case 2: // case 4
			pcase += (MC33_interiorTest(k, 0, v)? 239 + 6*k: 231 + (k<<1));
			break;
		case 3: // case 6
			if ((m? i: i^0xFF)&MC33_faceTest1(k%6, v))
				pcase += 575 + 5*k; //6.2
			else
				pcase += (MC33_interiorTest(k/6, 0, v)? 407 + 7*k: 335 + 3*k); //6.1
			break;
		case 4: // case 7
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -3:
				pcase += 695 + 3*k; //7.1
				break;
			case -1: //7.2
				pcase += (u.f[4] + u.f[5] < 0? (u.f[0] + u.f[2] < 0? 759: 799): 719) + 5*k;
				break;
			case 1: //7.3
				pcase += (u.f[4] + u.f[5] < 0? 983: (u.f[0] + u.f[2] < 0? 839: 911)) + 9*k;
				break;
			default: //7.4
				pcase += (MC33_interiorTest(k>>1, 0, v)? 1095 + 9*k: 1055 + 5*k);
			}
			break;
		case 5: // case 10
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2:
				if (k == 2? MC33_interiorTest(0, 0, v): MC33_interiorTest(0, 0, v)||MC33_interiorTest(k? 1: 3, 0, v))
					pcase += 1213 + (k<<3); //10.1.2
				else
					pcase += 1189 + (k<<2); //10.1.1
				break;
			case 0: //10.2
				pcase += (u.f[2 + k] < 0? 1261: 1285) + (k<<3);
				break;
			default:
				if (k == 2? MC33_interiorTest(1, 0, v): MC33_interiorTest(2, 0, v)||MC33_interiorTest(k? 3: 1, 0, v))
					pcase += 1237 + (k<<3); //10.1.2
				else
					pcase += 1201 + (k<<2); //10.1.1
			}
			break;
		case 6: // case 12
			switch (MC33_faceTests(u.f,(m? i: i^0xFF), v)) {
			case -2: //12.1
				pcase += (MC33_interiorTest((0xDA010C>>(k<<1))&3, 0, v)? 1453 + (k<<3): 1357 + (k<<2));
				break;
			case 0: //12.2
				pcase += (u.f[k>>1] < 0? 1645: 1741) + (k<<3);
				break;
			default: //12.1
				pcase += (MC33_interiorTest((0xA7B7E5>>(k<<1))&3, 0, v)? 1549 + (k<<3): 1405 + (k<<2));
			}
			break;
		default: // case 13
			switch (abs(MC33_faceTests(u.f, 165, v))) {
			case 0:
				k = ((u.f[1] < 0)<<1)|(u.f[5] < 0);
				if (u.f[0]*u.f[1] == u.f[5]) //13.4
					pcase += 2157 + 12*k;
				else {
					c = MC33_interiorTest(k, 1, v); // 13.5.1 if c == 0 else 13.5.2
					pcase += 2285 - 40*c + (c? 10: 6)*k;
				}
				break;
			case 2: //13.3
				pcase += 1917 + 10*((u.f[0] < 0? u.f[2] > 0: 12 + (u.f[2] < 0)) + (u.f[1] < 0? u.f[3] < 0: 6 + (u.f[3] > 0)));
				if (u.f[4] > 0)
					pcase += 30;
				break;
			case 4: //13.2
				k = 21 + 11*u.f[0] + 4*u.f[1] + 3*u.f[2] + 2*u.f[3] + u.f[4];
				if (k >> 4)
					k -= (k&32? 20: 10);
				pcase += 1845 + 3*k;
				break;
			default: //13.1
				pcase += 1839 + 2*u.f[0];
			}
	}
	while (i) {
		i = *(++pcase);
		for (k = 3; k;) {
			c = i&0x0F;
			i >>= 4;
			if (p[c] == FF) {
				switch (c) {
				case 0:
					if (z || x)
						p[0] = M->Dy[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[0] = p[3];
							else if (p[8] != FF)
								p[0] = p[8];
							else if (y && signbf(v[3]))
								p[0] = M->Lz[y][0];
							else if (y && signbf(v[4]))
								p[0] = M->Dx[y][0];
							else if (y? signbf(M->iso - M->F[0][y - 1][0]): 0)
								p[0] = M->Dy[y - 1][0];
							else
								p[0] = M->nV++;
						} else if (v[1] == 0) {
							if (p[1] != FF)
								p[0] = p[1];
							else if (p[9] != FF)
								p[0] = p[9];
							else
								p[0] = M->nV++;
						} else
							p[0] = M->nV++;

						M->Dy[y][0] = p[0];
					}
					break;
				case 1:
					if (x)
						p[1] = M->Lz[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[1] = p[0];
							else if (p[9] != FF)
								p[1] = p[9];
							else if (z && signbf(v[0]))
								p[1] = M->Dy[y][0];
							else if (z && signbf(v[5]))
								p[1] = M->Dx[y + 1][0];
							else if (z && y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][0]): 0)
								p[1] = M->Dy[y + 1][0];
							else if (z? signbf(M->iso - M->F[z - 1][y + 1][0]): 0) {
								p[1] = M->Lz[y + 1][0];
								break;
							} else
								p[1] = M->nV++;
						} else if (v[2] == 0) {
							if (p[2] != FF)
								p[1] = p[2];
							else if (p[10] != FF)
								p[1] = p[10];
							else
								p[1] = M->nV++;
						} else
							p[1] = M->nV++;
						M->Lz[y + 1][0] = p[1];
					}
					break;
				case 2:
					if (x)
						p[2] = M->Uy[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[2] = p[3];
							else if (p[11] != FF)
								p[2] = p[11];
							else if (y && signbf(v[0]))
								p[2] = M->Lz[y][0];
							else if (y && signbf(v[7]))
								p[2] = M->Ux[y][0];
							else if (y? signbf(M->iso - M->F[z + 1][y - 1][0]): 0)
								p[2] = M->Uy[y - 1][0];
							else
								p[2] = M->nV++;
						} else if (v[2] == 0) {
							if (p[1] != FF)
								p[2] = p[1];
							else if (p[10] != FF)
								p[2] = p[10];
							else
								p[2] = M->nV++;
						} else
							p[2] = M->nV++;
						M->Uy[y][0] = p[2];
					}
					break;
				case 3:
					if (y || x)
						p[3] = M->Lz[y][x];
					else {
						if (v[0] == 0) {
							if (p[0] != FF)
								p[3] = p[0];
							else if (p[8] != FF)
								p[3] = p[8];
							else if (z && signbf(v[1]))
								p[3] = M->Dy[0][0];
							else if (z && signbf(v[4]))
								p[3] = M->Dx[0][0];
							else if (z? signbf(M->iso - M->F[z - 1][0][0]): 0) {
								p[3] = M->Lz[0][0];
								break;
							} else
								p[3] = M->nV++;
						} else if (v[3] == 0) {
							if (p[2] != FF)
								p[3] = p[2];
							else if (p[11] != FF)
								p[3] = p[11];
							else
								p[3] = M->nV++;
						} else
							p[3] = M->nV++;
						M->Lz[0][0] = p[3];
					}
					break;
				case 4:
					if (z)
						p[4] = M->Dy[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[4] = p[8];
							//else if (p[7] != FF)
							//	p[4] = p[7];
							else if (y && signbf(v[0]))
								p[4] = M->Dx[y][x];
							else if (y && signbf(v[7]))
								p[4] = M->Lz[y][x + 1];
							else if (y? signbf(M->iso - M->F[0][y - 1][x + 1]): 0)
								p[4] = M->Dy[y - 1][x + 1];
							else if (y && x + 1 < M->nx? signbf(M->iso - M->F[0][y][x + 2]): 0)
								p[4] = M->Dx[y][x + 1];
							else
								p[4] = M->nV++;
						} else if (v[5] == 0) {
							if (p[9] != FF)
								p[4] = p[9];
							else if (p[5] != FF)
								p[4] = p[5];
							else
								p[4] = M->nV++;
						} else
							p[4] = M->nV++;
						M->Dy[y][x + 1] = p[4];
					}
					break;
				case 5:
					if (v[5] == 0) {
						if (z) {
							if (signbf(v[4]))
								p[5] = p[4] = M->Dy[y][x + 1];
							else if (signbf(v[1]))
								p[5] = p[9] = M->Dx[y + 1][x];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z][y + 1][x + 2]): 0)
								p[5] = M->Dx[y + 1][x + 1];
							else if (y + 1 < M->ny? signbf(M->iso - M->F[z][y + 2][x + 1]): 0)
								p[5] = M->Dy[y + 1][x + 1];
							else if (signbf(M->iso - M->F[z - 1][y + 1][x + 1])) {
								p[5] = M->Lz[y + 1][x + 1]; // value of previous slice
								break;
							} else
								p[5] = M->nV++;
						} else
							p[5] = M->nV++;
					} else
						p[5] = M->nV++;
					M->Lz[y + 1][x + 1] = p[5];
					break;
				case 6:
					if (v[7] == 0) {
						if (y) {
							if (signbf(v[3]))
								p[6] = p[11] = M->Ux[y][x];
							else if (signbf(v[4]))
								p[6] = p[7] = M->Lz[y][x + 1];
							else if (signbf(M->iso - M->F[z + 1][y - 1][x + 1]))
								p[6] = M->Uy[y - 1][x + 1];
							else if (x + 1 < M->nx? signbf(M->iso - M->F[z + 1][y][x + 2]): 0)
								p[6] = M->Ux[y][x + 1];
							else
								p[6] = M->nV++;
						} else if (p[11] != FF)
								p[6] = p[11];
							//else if (p[7] != FF)
							//	p[6] = p[7];
							else
								p[6] = M->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[6] = (p[10] == FF? M->nV++: p[10]);
						else
							p[6] = p[5];
					} else
						p[6] = M->nV++;
					M->Uy[y][x + 1] = p[6];
					break;
				case 7:
					if (y)
						p[7] = M->Lz[y][x + 1];
					else {
						if (v[4] == 0) {
							if (p[8] != FF)
								p[7] = p[8];
							else if (p[4] != FF)
								p[7] = p[4];
							else if (z && signbf(v[0]))
								p[7] = M->Dx[0][x];
							else if (z && signbf(v[5]))
								p[7] = M->Dy[0][x + 1];
							else if (z && x + 1 < M->nx? signbf(M->iso - M->F[z][0][x + 2]): 0)
								p[7] = M->Dx[0][x + 1];
							else if (z? signbf(M->iso - M->F[z - 1][0][x + 1]): 0) {
								p[7] = M->Lz[0][x + 1];
								break;
							} else
								p[7] = M->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[7] = p[6];
							else if (p[11] != FF)
								p[7] = p[11];
							else
								p[7] = M->nV++;
						} else
							p[7] = M->nV++;
						M->Lz[0][x + 1] = p[7];
					}
					break;
				case 8:
					if (z || y)
						p[8] = M->Dx[y][x];
					else {
						if (v[0] == 0) {
							if (p[3] != FF)
								p[8] = p[3];
							else if (p[0] != FF)
								p[8] = p[0];
							else if (x && signbf(v[3]))
								p[8] = M->Lz[0][x];
							else if (x && signbf(v[1]))
								p[8] = M->Dy[0][x];
							else if (x? signbf(M->iso - M->F[0][0][x - 1]): 0)
								p[8] = M->Dx[0][x - 1];
							else
								p[8] = M->nV++;
						} else if (v[4] == 0) {
							if (p[4] != FF)
								p[8] = p[4];
							else if (p[7] != FF)
								p[8] = p[7];
							else
								p[8] = M->nV++;
						} else
							p[8] = M->nV++;
						M->Dx[0][x] = p[8];
					}
					break;
				case 9:
					if (z)
						p[9] = M->Dx[y + 1][x];
					else {
						if (v[1] == 0) {
							if (p[0] != FF)
								p[9] = p[0];
							else if (p[1] != FF)
								p[9] = p[1];
							else if (x && signbf(v[0]))
								p[9] = M->Dy[y][x];
							else if (x && signbf(v[2]))
								p[9] = M->Lz[y + 1][x];
							else if (x? signbf(M->iso - M->F[0][y + 1][x - 1]): 0)
								p[9] = M->Dx[y + 1][x - 1];
							else
								p[9] = M->nV++;
						} else if (v[5] == 0) {
							if (p[5] != FF)
								p[9] = p[5];
							else if (p[4] != FF)
								p[9] = p[4];
							else
								p[9] = M->nV++;
						} else
							p[9] = M->nV++;
						M->Dx[y + 1][x] = p[9];
					}
					break;
				case 10:
					if (v[2] == 0) {
						if (x) {
							if (signbf(v[1]))
								p[10] = p[1] = M->Lz[y + 1][x];
							else if (signbf(v[3]))
								p[10] = p[2] = M->Uy[y][x];
							else if (signbf(M->iso - M->F[z + 1][y + 1][x - 1]))
								p[10] = M->Ux[y + 1][x - 1];
							else
								p[10] = M->nV++;
						} else if (p[2] != FF)
								p[10] = p[2];
							//else if (p[1] != FF)
							//	p[10] = p[1];
							else
								p[10] = M->nV++;
					} else if (v[6] == 0) {
						if (p[5] == FF)
							p[10] = (p[6] == FF? M->nV++: p[6]);
						else
							p[10] = p[5];
					} else
						p[10] = M->nV++;
					M->Ux[y + 1][x] = p[10];
					break;
				case 11:
					if (y)
						p[11] = M->Ux[y][x];
					else {
						if (v[3] == 0) {
							if (p[3] != FF)
								p[11] = p[3];
							else if (p[2] != FF)
								p[11] = p[2];
							else if (x && signbf(v[0]))
								p[11] = M->Lz[0][x];
							else if (x && signbf(v[2]))
								p[11] = M->Uy[0][x];
							else if (x? signbf(M->iso - M->F[z + 1][0][x - 1]): 0)
								p[11] = M->Ux[0][x - 1];
							else
								p[11] = M->nV++;
						} else if (v[7] == 0) {
							if (p[6] != FF)
								p[11] = p[6];
							else if (p[7] != FF)
								p[11] = p[7];
							else
								p[11] = M->nV++;
						} else
							p[11] = M->nV++;
						M->Ux[0][x] = p[11];
					}
				break;
				default:
					p[12] = M->nV++;
				}
			}
			u.ti[--k] = p[c];
		}
		if (u.ti[0] != u.ti[1] && u.ti[0] != u.ti[2] && u.ti[1] != u.ti[2])
			M->nT++;
	}
}
#undef FF

void MC33_freeTemp_O_N(MC33 *M) {
	free(M->Dx); free(M->Ux); free(M->Dy); free(M->Uy);
	free(M->Lz);
	free(M);
}

void free_MC33(MC33 *M) {
	unsigned int y;
	if (M) {
		for (y = 0; y != M->ny; y++) {
			free(M->Dx[y]); free(M->Ux[y]); free(M->Dy[y]); free(M->Uy[y]);
			free(M->Lz[y]);
		}
		free(M->Dx[M->ny]); free(M->Ux[M->ny]);
		free(M->Lz[M->ny]);
		MC33_freeTemp_O_N(M);
	}
}

MC33 *create_MC33(_GRD* G) {
	unsigned int x, y;
	MC33 *M;
	if (!G)
		return 0;
	M = (MC33*)malloc(sizeof(MC33));
	if (!M)
		return 0;
	M->nx = G->N[0];
	M->ny = G->N[1];
	M->nz = G->N[2];

#ifndef GRD_ORTHOGONAL
	if (G->nonortho) {
		M->store = MC33_spnC;
		for (int j = 0; j != 3; j++)
			for (int i = 0; i != 3; i++) {
				M->_A[j][i] = G->_A[j][i]*G->d[i]; // true transformation matrices
				M->A_[j][i] = G->A_[j][i]/G->d[j];
			}
	} else
#endif
	if (G->d[0] != G->d[1] || G->d[1] != G->d[2]) {
		M->ca = G->d[2]/G->d[0];
		M->cb = G->d[2]/G->d[1];
		M->store = MC33_spnB;
	} else
		M->store = (G->d[0] == 1 && G->r0[0] == 0 && G->r0[1] == 0 && G->r0[2] == 0? MC33_spn0: MC33_spnA);

	for (int j = 0; j != 3; j++) {
		M->O[j] = G->r0[j];
		M->D[j] = G->d[j];
	}
	M->Lz = (unsigned int**)malloc((M->ny + 1)*sizeof(int*));//edges 1, 3, 5 (only write) and 7
	M->Dy = (unsigned int**)malloc(M->ny*sizeof(int*));//edges 0 (only read) and 4
	M->Uy = (unsigned int**)malloc(M->ny*sizeof(int*));//edges 2 and 6 (only write)
	M->Dx = (unsigned int**)malloc((M->ny + 1)*sizeof(int*));//edges 8 and 9
	M->Ux = (unsigned int**)malloc((M->ny + 1)*sizeof(int*));//edges 10 (only write) and 11
	if (!M->Ux) {
		MC33_freeTemp_O_N(M);
		return 0;
	}
	M->F = (const GRD_data_type***)G->F;
	x = M->nx*sizeof(int);
	for (y = 0; y != M->ny; y++) {
		M->Dx[y] = (unsigned int*)malloc(x);
		M->Ux[y] = (unsigned int*)malloc(x);
		M->Lz[y] = (unsigned int*)malloc(x + sizeof(int));
		M->Dy[y] = (unsigned int*)malloc(x + sizeof(int));
		M->Uy[y] = (unsigned int*)malloc(x + sizeof(int));
	}
	if (M->Uy[y - 1]) {
		M->Dx[y] = (unsigned int*)malloc(x);
		M->Ux[y] = (unsigned int*)malloc(x);
		M->Lz[y] = (unsigned int*)malloc(x + sizeof(int));
		if (M->Lz[y])
			return M;
	} else
		M->Dx[y] = M->Ux[y] = M->Lz[y] = 0;
	free_MC33(M);
	return 0;
}

surface* calculate_isosurface(MC33 *M, float iso) {
	unsigned int x, y, z, Nx = M->nx;
	float Vt[12];
	float *v1 = Vt, *v2 = Vt + 4;
	const GRD_data_type ***F = M->F, **F0, **F1, *V00, *V01, *V11, *V10;
	surface *S = (surface*)malloc(sizeof(surface));
	if (!S)
		return 0;
	M->nT = M->nV = 0;
	M->memoryfault = 0;
	M->capt = M->capv = 4096;
	M->T = (unsigned int(*)[3])malloc(3*4096*sizeof(int));
	M->N = (float(*)[3])malloc(3*4096*sizeof(float));
	M->V = (float(*)[3])malloc(3*4096*sizeof(float));
	M->iso = iso;
	if (M->V)
		for (z = 0; z != M->nz; z++) {
			F0 = *F;
			F1 = *(++F);
			for (y = 0; y != M->ny; y++) {
				V00 = *F0;
				V01 = *(++F0);
				V10 = *F1;
				V11 = *(++F1);
				v2[0] = iso - *V00;//the difference was inverted to use signbit function
				v2[1] = iso - *V01;
				v2[2] = iso - *V11;
				v2[3] = iso - *V10;
				//the eight least significant bits of i correspond to vertex indices. (x...x01234567)
				//If the bit is 1 then the vertex value is greater than zero.
				unsigned int i = signbf(v2[3]) != 0;
				if (signbf(v2[2])) i |= 2;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[0])) i |= 8;
				for (x = 0; x != Nx; x++) {
					{float *P = v1; v1 = v2; v2 = P;}//v1 and v2 are exchanged
					v2[0] = iso - *(++V00);
					v2[1] = iso - *(++V01);
					v2[2] = iso - *(++V11);
					v2[3] = iso - *(++V10);
					i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
					if (signbf(v2[2])) i |= 2;
					if (signbf(v2[1])) i |= 4;
					if (signbf(v2[0])) i |= 8;
					if (i && i^0xFF) {
						if (v1 > v2) {float *t = v2; float *s = t + 8; *s = *t; *(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t);}
						MC33_findCase(M,x,y,z,i,v1);
					}
				}
			}
			{unsigned int** P = M->Dx; M->Dx = M->Ux; M->Ux = P;}//M->Dx and M->Ux are exchanged
			{unsigned int** P = M->Dy; M->Dy = M->Uy; M->Uy = P;}//M->Dy and M->Uy are exchanged
		}
	else
		M->memoryfault = 1;
	if (M->nV) {
		M->color = (int*)malloc(M->nV*sizeof(int));
		memcpy(S, M, offsetof(MC33, iso));
		if (M->color) {
			*M->color = DefaultColorMC;
			while (--M->nV)
				*(++M->color) = DefaultColorMC;
		} else
			M->memoryfault = 1;
	} else {
		free(M->V); free(M->N); free(M->T);
		memset(S, 0, sizeof(surface));
	}
	if (M->memoryfault) {
		free_surface_memory(S);
		return 0;
	}
	S->iso = iso;
	return S;
}

// modified from calculate_isosurface function
unsigned long long size_of_isosurface(MC33 *M, float iso, unsigned int *nV, unsigned int *nT) {
	unsigned int x, y, z, Nx = M->nx;
	float Vt[12];
	float *v1 = Vt, *v2 = Vt + 4;
	const GRD_data_type ***F = M->F, **F0, **F1, *V00, *V01, *V11, *V10;
	M->nT = M->nV = 0;
	M->iso = iso;
	for (z = 0; z != M->nz; z++) {
		F0 = *F;
		F1 = *(++F);
		for (y = 0; y != M->ny; y++) {
			V00 = *F0;
			V01 = *(++F0);
			V10 = *F1;
			V11 = *(++F1);
			v2[0] = iso - *V00;
			v2[1] = iso - *V01;
			v2[2] = iso - *V11;
			v2[3] = iso - *V10;
			unsigned int i = signbf(v2[3]) != 0;
			if (signbf(v2[2])) i |= 2;
			if (signbf(v2[1])) i |= 4;
			if (signbf(v2[0])) i |= 8;
			for (x = 0; x != Nx; x++) {
				{float *P = v1; v1 = v2; v2 = P;}
				v2[0] = iso - *(++V00);
				v2[1] = iso - *(++V01);
				v2[2] = iso - *(++V11);
				v2[3] = iso - *(++V10);
				i = ((i&0x0F)<<4)|(signbf(v2[3]) != 0);
				if (signbf(v2[2])) i |= 2;
				if (signbf(v2[1])) i |= 4;
				if (signbf(v2[0])) i |= 8;
				if (i && i^0xFF) {
					if (v1 > v2) {float *t = v2; float *s = t + 8; *s = *t; *(++s) = *(++t); *(++s) = *(++t); *(++s) = *(++t);}
					MC33_Case_count(M,x,y,z,i,v1);
				}
			}
		}
		{unsigned int** P = M->Dx; M->Dx = M->Ux; M->Ux = P;}
		{unsigned int** P = M->Dy; M->Dy = M->Uy; M->Uy = P;}
	}
	if (nV)
		*nV = M->nV;
	if (nT)
		*nT = M->nT;
	// number of vertices * (size of vertex and normal + size of color ) + number of triangle * size of triangle + size of struct surface
	return M->nV * (6 * sizeof(float) + sizeof(int)) + M->nT * (3 * sizeof(int)) + sizeof(surface);
}

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif //marching_cubes_33_c

