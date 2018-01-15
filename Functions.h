/* This file lists all the functions */

#include"GlobalVariables.h"
// read in the geometry and build the list.
// input: the input file name/address.
// output: None.
int GeoBuilder(char *fname);

// free all the allocated memory.
// input: argc.
// output: none.
void FreeMem(int argc);

// calculate the thermal conductivity.
// input:The bulk MFP spectra.
// output: effective thermal conductivity to file.
int CalTherm(char *fname);

// calculate the thermal conductivity.
// input:The bulk MFP spectra.
// output: effective thermal conductivity to file.
void CalMFPs(char *fname);

// Calculate 1~20 Periods both specularity of 0 and 1.
// input:filename
// output: MFP_bdy
void CalMFPs(char *pathgeo);

// Calculate the transmission
// input: bulkMFP
// output: the transmission rate
double CalTransmission(double bulkMFP);

// Perform Ray Tracing.
// input: bulkMFP, particle
// output: reflected or transmitted?
void PerformRayTracing(double bulkMFP, particle p);

// build cumulative function for choosing starting surface
// input: none; prerequisite: geobuilder;
// output: none; Area CF is computed
void buildCumuAreas(void);

// tell if the point is inside an triangle in 2D
// input: Surface Number, point position.
// output: true or false: 1 or 0;
int IsInsideTriangle2D(double px, double py, int SurfaceNum);

// Random number Generator
// input: none
// output: a uniformly distributed random number in [0,1]
double RandR(void);

// RNG Initializer
// input: seed
// output: none
void InitRand(unsigned long long Seed);

// Calculate the area of triangle.
// input: three vertices
// output: area
double TriArea(double p1[3], double p2[3], double p3[3]);

// Other claimers:
// The default unit in this code is nm for convenience.
