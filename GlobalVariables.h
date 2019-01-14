/* This file stores the global variables */
// Avoid redefine
#ifndef VAR_H
#define VAR_H
// Pointers:
extern double **Points;                          // N x 3 array;
extern int **Surfaces;                           // N x 3 array; integer point to the vertices.
extern double **SurfaceNorms;                    // N x 3 array;
extern double *SurfaceAreas;                     // N x 1 array;
extern double *Specularity;                      // N x 1 array;
extern int *InSurfaces;                          // k x 1 array; integer point to the surfaces.
extern int *OutSurfaces;                         // k x 1 array; integer point to the surfaces.
extern double *Displacement;                     // N x 1 array; for calculating distance.
extern int *SurfaceMarker;                       // N x 1 array; mark the surface types.
extern double *CumuAreas;                        // N x 1 array; Surface Cumulative Areas for choosing the starting surface.
extern double *dist2surf;                        // N x 1 array: distance to surface.

// MFP Spectra
extern int NumBins;                              // recording the MFP bin number
extern double *MFPs;                             // recording the MFPs
extern double *Klambda;                          // recording the Klambda;
extern double *dlambda;                          // recording the bin width.

// Final results
extern double effectiveTC;                       // Effective Thermal Conductivity.    Need to intialize once.

// Rand Seed
extern unsigned long long RandSeed;              // Random number generator seed;

// Calculation Parameters:
extern int NumPoints;                            // Number of Points;
extern int NumSurfaces;                          // Number of Surfaces;
extern int NumPeriods;                           // Simulation period number;
extern int NumOut;                               // Number of Transmitted Particles.
extern int NumInSurfaces;                        // Number of Surfaces in inflow.
extern int NumOutSurfaces;                       // Number of Surfaces in outflow.
extern int NumParticles;                         // Number of Particles in simulation.re-intialize for different bulkMFP
extern int NumTrans;                             // Number of transmitted Particles.  re-intialize for different bulkMFP
extern int NumReflect;                           // Number of reflected Particles.    re-intialize for different bulkMFP
extern double UnitCellLength;                      // Input the sample length.Thus the inflow and outflow surfaces are parallel.
extern double Spec;                              // Input the specularity of the outer surfaces;

// Structure;
typedef struct particle{
	double x;                             // location @ x;
	double y;                             // location @ y;
	double z;                             // location @ z;
	double tx;                            // direction @ x;
	double ty;                            // direction @ y;
	double tz;                            // direction @ z;
	int Surface;                          // collision surface.
}particle;

// Fixed Constant;
extern const double INF;                       // Represent the infinity.
extern const double PI;                        // Constant PI.
#endif
