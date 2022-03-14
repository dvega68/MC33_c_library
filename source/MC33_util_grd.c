/*
	File: MC33_util_grd.c
	Programed by: David Vega: dvega@uc.edu.ve
	March 2012
	updated by David Vega
	November 2017.
	February 2018.
	April 2019.
	July 2019.
	February 2020.
	August 2020.
	August 2021.
	December 2021.
*/

#ifndef marching_cubes_33_h
#error ***Include the file 'marching_cubes_33.h' or 'marching_cubes_33.c' before***
#include "*" //to abort the compilation
#endif

#ifndef MC33_util_grd_c
#define MC33_util_grd_c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef GRD_orthogonal
void setIdentMat3x3d(double (*A)[3]) {
	for (double *d = A[0] + 8; --d != A[0];)
		d[0] = 0.0;
	for (int i = 0; i != 3; ++i)
		A[i][i] = 1.0;
}
#endif

#ifndef LAGRANGE3D4GRD_H
//******************************************************************
void free_memory_grd(_GRD *Z) {
	unsigned int k, j;
	if (Z) {
		if (Z->F) {
			if (Z->internal_data)
				for (k = 0; k <= Z->N[2]; ++k) {
					if (Z->F[k]) {
						for (j = 0; j <= Z->N[1]; ++j)
							free(Z->F[k][j]);
					} else {
						k = Z->N[2];
						break;
					}
					free(Z->F[k]);
				}
			else
				for (k = 0; k <= Z->N[2]; ++k)
					free(Z->F[k]);
			free(Z->F);
		}
		free(Z);
	}
}

int alloc_F(_GRD* Z) {
	unsigned int j, k;
	Z->F = (GRD_data_type***)malloc((Z->N[2] + 1)*sizeof(void*));
	if (!Z->F)
		return -1;
	for (k = 0; k <= Z->N[2]; ++k) {
		Z->F[k] = (GRD_data_type**)malloc((Z->N[1] + 1)*sizeof(void*));
		if (!Z->F[k])
			return -1;
		for (j = 0; j <= Z->N[1]; ++j) {
			Z->F[k][j] = (GRD_data_type*)malloc((Z->N[0] + 1)*sizeof(GRD_data_type));
			if (!Z->F[k][j]) {
				while (j)
					free(Z->F[k][--j]);
				free(Z->F[k]);
				Z->F[k] = 0;
				return -1;
			}
		}
	}
	Z->internal_data = 1;
	return 0;
}


//******************************************************************
/*
read_grd read a filename file (the file must be a output *.grd file from the
DMol program), it returns a pointer to struct _GRD that contains all the grid
data.
*/
_GRD* read_grd(const char *filename) {
	_GRD *Z;
	FILE *in;
	char line[128];
	unsigned int i, j, k;
	double ca, cb, sg, cg, aux1, aux2;
	int grd_xi[3], grd_ordenij;
	float Ang[3];

	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;

#ifndef GRD_orthogonal
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	in = fopen(filename,"r");
	if (!in) return 0;
	fgets(Z->title,159,in);
	fgets(line,60,in);
	fgets(line,60,in);
	sscanf(line,"%f %f %f %f %f %f",&(Z->L[0]),&(Z->L[1]),&(Z->L[2]),&Ang[0],&Ang[1],&Ang[2]);
	fgets(line,60,in);
	sscanf(line,"%d %d %d",&(Z->N[0]),&(Z->N[1]),&(Z->N[2]));
	fgets(line,60,in);
	sscanf(line,"%d %d %*d %d %*d %d %*d",&grd_ordenij,&grd_xi[0],&grd_xi[1],&grd_xi[2]);
	if (Z->N[0] < 2 || Z->N[1] < 2 || Z->N[2] < 2) return 0;
	if (grd_ordenij != 1 && grd_ordenij != 3) return 0;
	for (i = 0; i != 3; ++i) {
		Z->d[i] = Z->L[i]/Z->N[i];
		Z->r0[i] = grd_xi[i]*Z->d[i];
	}

	Z->periodic = (grd_xi[0] == 0)|((grd_xi[1] == 0)<<1)|((grd_xi[2] == 0)<<2);

#ifndef GRD_orthogonal
	memcpy(Z->Ang, Ang, sizeof Ang);
	if (Ang[0] != 90 || Ang[1] != 90 || Ang[2] != 90) {
		Z->nonortho = 1;
		ca = cos(Ang[0]*(M_PI/180.0));
		cb = cos(Ang[1]*(M_PI/180.0));
		aux1 = Ang[2]*(M_PI/180.0);
		sg = sin(aux1);
		cg = cos(aux1);
		aux1 = ca - cb*cg;
		aux2 = sqrt(sg*sg + 2*ca*cb*cg - ca*ca - cb*cb);
		Z->_A[0][0] = Z->A_[0][0] = 1.0;
		Z->_A[0][1] = cg;
		Z->_A[0][2] = cb;
		Z->_A[1][1] = sg;
		Z->A_[1][1] = cb = 1.0/sg;
		Z->A_[0][1] = -cg*cb;
		Z->_A[1][2] = aux1*cb;
		Z->_A[2][2] = aux2*cb;
		aux2 = 1.0/aux2;
		Z->A_[0][2] = (cg*aux1 - ca*sg*sg)*cb*aux2;
		Z->A_[1][2] = -aux1*cb*aux2;
		Z->A_[2][2] = sg*aux2;
	} else {
		Z->nonortho = 0;
		setIdentMat3x3d(Z->A_);
		setIdentMat3x3d(Z->_A);
	}
#endif

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

	for (k = 0; k <= Z->N[2]; ++k)
		if (grd_ordenij == 1) {
			for (j = 0; j <= Z->N[1]; ++j)
				for (i = 0; i <= Z->N[0]; ++i)
					fscanf(in,"%f",&Z->F[k][j][i]);
		} else {
			for (i = 0; i <= Z->N[0]; ++i)
				for (j = 0; j <= Z->N[1]; ++j)
					fscanf(in,"%f",&Z->F[k][j][i]);
		}
	fclose(in);
	return Z;
}

/*
internal binary format
*/
_GRD* read_grd_binary(const char* filename) {
	FILE* in;
	unsigned int i, j, k;
	_GRD* Z;

	in = fopen(filename,"rb");
	if (!in)
		return 0;
	fread(&i,sizeof(int),1,in);
	if (i != 0x4452475f) // _GRD
		return 0;

	fread(&i,sizeof(int),1,in);
	if (i > 159)
		return 0;
	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;
	Z->internal_data = 1;
	fread(Z->title,i*sizeof(char),1,in);
//	Z->title[i] = '\0';
	fread(Z->N,sizeof Z->N,1,in);
	fread(Z->L,sizeof Z->L,1,in);
	fread(Z->r0,sizeof Z->r0,1,in);
	fread(Z->d,sizeof Z->d,1,in);
#ifndef GRD_orthogonal
	fread(&Z->nonortho,sizeof(int),1,in);
	if (Z->nonortho) {
		fread(Z->Ang,3*sizeof(float),1,in);
		fread(Z->_A,sizeof Z->_A,1,in);
		fread(Z->A_,sizeof Z->A_,1,in);
		mult_Abf = _multA_bf;
	} else {
		setIdentMat3x3d(Z->A_);
		setIdentMat3x3d(Z->_A);
	}
#else
	fread(&i, sizeof(int), 1, in);
	if (i)
		fseek(in, 3*sizeof(float) + 18*sizeof(double), SEEK_CUR);
#endif
	//Z->periodic = 0;
	if (Z->r0[0] == 0 && Z->r0[1] == 0 && Z->r0[2] == 0)
		Z->periodic = 1;

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

	for (k = 0; k <= Z->N[2]; ++k)
		for (j = 0; j <= Z->N[1]; ++j)
			fread(Z->F[k][j],(Z->N[0] + 1)*sizeof(GRD_data_type),1,in);
	fclose(in);
	return Z;
}
#endif
//******************************************************************
/*
Reads a set of files that contain a slab of res*res scan data points, the data
points are read as unsigned short int (if order is different from 0, the bytes
of the unsigned short are exchanged). The filename must end with a number, and
the fuction read all files with end number greater or equal to filename.
(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)
*/
_GRD* read_scanfiles(const char *filename, unsigned int res, int order) {
	_GRD *Z;
	FILE *in;
	char *nm, *nm2;
	unsigned int i, j, l;
	int k = -1;
	unsigned short int n;
	GRD_data_type ***pt, **p;

	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;

#ifndef GRD_orthogonal
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	memset(Z->r0,0,sizeof Z->r0);
	for (i = 0; i != 2; ++i) {
		Z->d[i] = 1;
		Z->L[i] = Z->N[i] = res - 1;
	}
	Z->d[2] = 1;
	Z->F = 0;
	l = strlen(filename) - 1;
	nm = (char*)malloc((l + 5)*sizeof(char));
	if (!nm)
		return 0;
	while (filename[l] >= '0' && filename[l] <= '9')
		l--;
	strncpy(nm,filename,l + 1);
	nm2 = nm + l + 1;
	l = atoi(filename + l + 1);
	while (1) {
		sprintf(nm2,"%-d",l++);
		in = fopen(nm,"rb");
		if (!in)
			break;
		if (!((++k)&63)) {
			pt = (GRD_data_type***)realloc(Z->F,(k + 64)*sizeof(void*));
			if(!pt)
				break;
			Z->F = pt;
		}
		Z->F[k] = (GRD_data_type**)malloc(res*sizeof(void*));
		if (!Z->F[k])
			break;
		for (j = 0; j != res; ++j)
			Z->F[k][j] = (GRD_data_type*)malloc(res*sizeof(GRD_data_type));
		if (!Z->F[k][Z->N[1]])
			break;

		if (order)
			for (j = 0; j != res; ++j)
				for (i = 0; i != res; ++i) {
					fread(&n,sizeof(short int),1,in);
					Z->F[k][j][i] = (unsigned short int)((n>>8)|(n<<8));
				}
		else
			for (j = 0; j != res; ++j)
				for (i = 0; i != res; ++i) {
					fread(&n,sizeof(short int),1,in);
					Z->F[k][j][i] = n;
				}
		fclose(in);
	}
	free(nm);
	Z->L[2] = Z->N[2] = k;
	j = k>>1;
	for (i = 0; i != j; ++i) {
		p = Z->F[i];
		Z->F[i] = Z->F[k - i];
		Z->F[k - i] = p;
	}
#ifndef GRD_orthogonal
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

//******************************************************************
/*
Reads a file that contains only the data points as integers of 8, 16 or 32 bits.
byte is the number of bytes of each data point (1 to 4). If the data is big endian,
byte must be negative. The vector N[3] contains the number of points in each
dimension. The size of file must be byte*N[0]*N[1]*N[2].
*/
_GRD* read_raw_file(const char *filename, unsigned int *N, int byte, int isfloat) {
	unsigned int i, j, k;
	_GRD *Z;
	FILE *in;
	unsigned int ui = 0;
	if (isfloat) {
		if (byte != 4 && byte != 8)
			return 0;
	} else if (abs(byte) > 4 || abs(byte) == 3 || !byte)
		return 0;
	if (byte == -1) byte = 1;
	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z)
		return 0;

#ifndef GRD_orthogonal
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	memset(Z->r0,0,sizeof Z->r0);
	for (i = 0; i != 3; ++i)
		Z->d[i] = 1;
	in = fopen(filename,"rb");
	if (!in)
		return 0;
	Z->L[0] = Z->N[0] = N[0] - 1;
	Z->L[1] = Z->N[1] = N[1] - 1;
	Z->L[2] = Z->N[2] = N[2] - 1;

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

#if defined(integer_GRD)
	if (!isfloat && size_type_GRD == byte)
#else
	if (isfloat && size_type_GRD == byte)
#endif
	{
		byte *= N[0];
		for (k = 0; k != N[2]; ++k)
			for (j = 0; j != N[1]; ++j)
				fread(Z->F[k][j],byte,1,in);
	} else if (isfloat) {
		if (byte == 8) {
#if defined(integer_GRD) || size_type_GRD == 4
			double df;
			for (k = 0; k != N[2]; ++k)
				for (j = 0; j != N[1]; ++j)
					for (i = 0; i != N[0]; ++i) {
						fread(&df,byte,1,in);
						Z->F[k][j][i] = (GRD_data_type)df;
					}
#endif
		} else {
#if defined(integer_GRD) || size_type_GRD == 8
			float f;
			for (k = 0; k != N[2]; ++k)
				for (j = 0; j != N[1]; ++j)
					for (i = 0; i != N[0]; ++i) {
						fread(&f,byte,1,in);
						Z->F[k][j][i] = (GRD_data_type)f;
					}
#endif
		}
	} else if (byte < 0) {
		byte = -byte;
		for (k = 0; k != N[2]; ++k)
		//for (k = N[2] - 1; k >= 0; --k)
			for (j = 0; j != N[1]; ++j)
				for (i = 0; i != N[0]; ++i) {
					fread(&ui,byte,1,in);
					if (byte == 2)
						Z->F[k][j][i] = (ui>>8)|((ui<<8)&0xff00);
					else //if (byte == 4)
						Z->F[k][j][i] = (ui>>24)|((ui>>8)&0xff00)|((ui<<8)&0xff0000)|(ui<<24);
				}
	} else {
		for (k = 0; k != N[2]; ++k)
		//for (k = N[2] - 1; k >= 0; --k)
			for (j = 0; j != N[1]; ++j)
				for (i = 0; i != N[0]; ++i) {
					fread(&ui,byte,1,in);
					Z->F[k][j][i] = ui;
				}
	}
	fclose(in);
#ifndef GRD_orthogonal
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

//******************************************************************
/*
Reads a dat file:
http://www.cg.tuwien.ac.at/research/vis/datasets/
*/
_GRD* read_dat_file(const char *filename) {
	_GRD *Z;
	FILE *in;
	unsigned int i, j;
	unsigned short int n, nx, ny, nz;

	Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z)
		return 0;

#ifndef GRD_orthogonal
	Z->nonortho = 0;
#endif
	Z->periodic = 0;
	Z->internal_data = 1;
	memset(Z->r0,0,sizeof Z->r0);
	for (i = 0; i != 3; ++i)
		Z->d[i] = 1;
	in = fopen(filename,"rb");
	if (!in)
		return 0;
	fread(&nx,sizeof(short int),1,in);
	fread(&ny,sizeof(short int),1,in);
	fread(&nz,sizeof(short int),1,in);
	Z->L[0] = Z->N[0] = nx - 1;
	Z->L[1] = Z->N[1] = ny - 1;
	Z->L[2] = Z->N[2] = nz - 1;

	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}

	while (nz--)
		for (j = 0; j != ny; ++j)
			for (i = 0; i != nx; ++i) {
				fread(&n,sizeof(short int),1,in);
				Z->F[nz][j][i] = n;
			}
	fclose(in);
#ifndef GRD_orthogonal
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}

//******************************************************************
/*
	set_data_pointer creates a _GRD struct from an external data array. data
	must be stored with the nested inner loop running from i = 0 to Nx - 1
	and the outer loop from k = 0 to Nz - 1. The data will not be erased by
	free_memory_grd function. The function returns a pointer to the created
	struct.
*/
_GRD* grid_from_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data) {
	if (!data || Nx == 0 || Ny == 0 || Nz == 0)
		return 0;
	unsigned int i, j, k;
	_GRD *Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;
	Z->internal_data = 0;
	Z->F = (GRD_data_type***)malloc(Nz*sizeof(void*));
	if (!Z->F) {
		free(Z);
		return 0;
	}
	Z->N[0] = Nx - 1;
	Z->N[1] = Ny - 1;
	Z->N[2] = Nz - 1;
	for (k = 0; k < Nz; ++k) {
		Z->F[k] =(GRD_data_type**)malloc(Ny*sizeof(void*));
		if (!Z->F[k]) {
			while (k)
				free(Z->F[--k]);
			free(Z->F);
			free(Z);
			return 0;
		}
		for (j = 0; j < Ny; ++j)
			Z->F[k][j] = data + j*Nx;
		data += Ny*Nx;
	}
	for (i = 0; i != 3; ++i) {
		Z->L[i] = Z->N[i];
		Z->d[i] = 1.0;
		Z->r0[i] = 0.0;
#ifndef GRD_orthogonal
		Z->Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	Z->nonortho = 0;
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}


_GRD* generate_grid_from_fn(double xi, double yi, double zi, double xf, double yf, double zf, double dx, double dy, double dz, double (*fn)(double x, double y, double z)) {
	if (dx <= 0 || dy <= 0 || dz <= 0 || xi == xf || yi == yf || zi == zf)
		return 0;

	if (xi > xf) {
		double t = xi; xi = xf; xf = t;}
	if (xf - xi < dx)
		dx = xf - xi;
	if (yi > yf) {
		double t = yi; yi = yf; yf = t;}
	if (yf - yi < dy)
		dy = yf - yi;
	if (zi > zf) {
		double t = zi; zi = zf; zf = t;}
	if (zf - zi < dz)
		dz = zf - zi;
	_GRD *Z = (_GRD*)malloc(sizeof(_GRD));
	if (!Z) return 0;
	Z->N[0] = (int)((xf - xi)/dx + 0.5);
	Z->N[1] = (int)((yf - yi)/dy + 0.5);
	Z->N[2] = (int)((zf - zi)/dz + 0.5);
	if (alloc_F(Z)) {
		free_memory_grd(Z);
		return 0;
	}
	Z->d[0] = dx; Z->d[1] = dy; Z->d[2] = dz;
	Z->r0[0] = xi; Z->r0[1] = yi; Z->r0[2] = zi;
	if (fn) {
		unsigned int i, j, k;
		double x, y, z = zi;
		for (k = 0; k <= Z->N[2]; ++k) {
			y = yi;
			for (j = 0; j <= Z->N[1]; ++j) {
				x = xi;
				for (i = 0; i <= Z->N[0]; ++i) {
					Z->F[k][j][i] = (GRD_data_type)fn(x,y,z);
					x += dx;
				}
				y += dy;
			}
			z += dz;
		}
	}
	for (int i = 0; i != 3; ++i) {
		Z->L[i] = Z->N[i]*Z->d[i];
#ifndef GRD_orthogonal
		Z->Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	Z->nonortho = 0;
	setIdentMat3x3d(Z->_A);
	setIdentMat3x3d(Z->A_);
#endif
	return Z;
}


//******************************************************************
#endif //MC33_util_grd_c

