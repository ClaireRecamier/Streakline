#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

int stream(double *x0, double *v0, double *xm1, double *xm2, double *xm3, double *xp1, double *xp2, double *xp3, double *vm1, double *vm2, double *vm3, double *vp1, double *vp2, double *vp3, double *par, double *offset, int potential, int integrator, int N, int M, double mcli, double mclf, double rcl, double dt_, FILE *fpt, FILE *fpt1, FILE *fpt2);

int orbit(double *x0, double *v0, double *x1, double *x2, double *x3, double *v1, double *v2, double *v3, double *par, int potential, int integrator, int N, double dt_, double direction);
int back(double *x0, double *v0, double *xf, double *vf, double *par, double hp, int potential, int integrator, int N);

void dostep(double *x, double *v, double *par, int potential, double deltat, double sign);
void dostep1(double *x, double *v, double *par, int potential, double deltat, double sign);

void dostep_stream(double *xc, double *x, double *v, double *par, int potential, double Mcl, double deltat, double sign);
void dostep_rk(double *x, double *v, double *par, int potential, double deltat, double sign);

void force(double *x, double *a, double *par, int potential);
void force_plummer(double *x, double *a, double Mcl);
void initpar(int potential, double *par, double *apar);
double jacobi(double *x, double *v, double *par, int potential, double Mcl);

double len(double *x);
void n2t(double *x, double *x1, double *x2, double *x3, int i);
void t2n(double *x, double *x1, double *x2, double *x3, int i);
void t2t(double *x, double *y);

double ran1(long *idum);
double gasdev(long *idum);

// Simulation parameters
extern double Mcli, Mclf, Rcl, dt;

// Physical constants
#define G 6.67e-11		// Gravitational constant
#define yr 3.15569e7		// yr in seconds
#define kpc 3.08567758e19	// kpc in m
#define kpcau 2.063e+8
#define kms 1e3			// km/s in m/s
#define Msun 1.99e30		// Mass of the Sun
#define pi 3.141592653589793

// NR random generator definitions
#define Nnr 90
#define MAX 200
#define NR 1000
#define ts 5776.0
#define zxs 0.023
#define dnus 135.
#define nums 3050.

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
