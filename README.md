**RayTracing_ThermalConductivity**
`var x - 7;`
This code is used to calculate the the effective MFP for nanostructure.

Hard Code Part:
1. The inlet and outlet surface are parallel and they are on a plane which is normal to Z direction.
2. The geometry consists of triangular surfaces.

Example Input: remove all the lines starting with '#'
File1: for the MFP spectra
/# use "#" to comment
/# First line: total number of MFP bins
100
/# Follows: central MFP | K_Lambda | d_Lambda
5     100    10
15    100    10
25    100    10
35    100    10
...
File2: for the geometry
/# The very beginning: Random Seed
9443413163
/# specularity for outer surfaces
0.5
/# Points: Total Num.
/# Point Num + (X,Y,Z)
/# the order does not necessary to be kept.
100         
0   x0  y0  z0
1   x1  y1  z1
2   x2  y2  z2
3   x3  y3  z3
...
99  x99 y99 z99     

/# Surfaces:
/# Total Num.
/# Surface Num + (p1,p2,p3,specularity,surfacetypemarker(in/out/none))
/# The order does not necessary to be kept.
20              
0   p1  p2  p3  marker
1   p1  p2  p3  marker
2   p1  p2  p3  marker
3   p1  p2  p3  marker
...
19  p1  p2  p3  marker        

/# assign starting surfaces
/# total num
/# component
3
4
7
9

/# assign outflow surfaces
/# total num
/# component
0
3
4

/# give periodicity of the structure.
3
/# Give the sample length of the structure (nm)
200

/# total simulation particle.
400000

/# Again, remove all the lines with '#' for input.
The geometry input can be generated with COMSOL with some modification.


The codes can be modified freely.

Reference Methodology: F. Yang and C. Dames, Physical Review B 87 (2013).

my page: [elancyu](https://github.com/elancyu/RayTracing_ThermalConductivity)
