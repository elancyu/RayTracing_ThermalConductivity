/* This file stores the global variables */

// Pointers:
double **Points;                          // N x 3 array;
int **Surfaces;                           // N x 3 array; integer point to the vertices.
double **SurfaceNorms;                    // N x 3 array;
double *SurfaceAreas;                     // N x 1 array;
double *Specularity;                      // N x 1 array;
int *InSurfaces;                          // k x 1 array; integer point to the surfaces.
int *OutSurfaces;                         // k x 1 array; integer point to the surfaces.
double *Displacement;                     // N x 1 array; for calculating distance.
int *SurfaceMarker;                       // N x 1 array; mark the surface types.
double *CumuAreas;                        // N x 1 array; Surface Cumulative Areas for choosing the starting surface.
double *dist2surf;

// MFP Spectra
int NumBins;                              // recording the MFP bin number
double *MFPs;                             // recording the MFPs
double *Klambda;                          // recording the Klambda;
double *dlambda;                          // recording the bin width.

// Final results
double effectiveTC;                       // Effective Thermal Conductivity.    Need to intialize once.

// RandSeed
unsigned long long RandSeed;              // Random Seed

// Calculation Parameters:
int NumPoints;                            // Number of Points;
int NumSurfaces;                          // Number of Surfaces;
int NumPeriods;                           // Simulation period number;
int NumInSurfaces;                        // Number of Surfaces in inflow.
int NumOutSurfaces;                       // Number of Surfaces in outflow.
int NumParticles=10;                  // Number of Particles in simulation.re-intialize for different bulkMFP
int NumTrans;                             // Number of transmitted Particles.  re-intialize for different bulkMFP
int NumReflect;                           // Number of reflected Particles.    re-intialize for different bulkMFP
double UnitCellLength;                      // Input the sample length.Thus the inflow and outflow surfaces are parallel.
double Spec;                              // Input specularity.

// Fixed Constant;
double const INF=1e20;                         // Represent the infinity.
double const PI=3.141592653;                    // Constant PI from UIUC: truncated at 9th digit.
