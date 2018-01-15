/* This C will read in geometry then store it */
// input: geometry.in
// output: none.

/* description for input file: the surface is consist of points
and points are in triplet of real number. First block is points
Second block is surface. and Final block is periodicity.     */



#include"GlobalVariables.h"
#include"Functions.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

int GeoBuilder(char* fname)
{
	FILE *fpgeo;                // geometry input file.
	FILE *fpchecker;           // geometry checker.
	int i, j;                   // loop integer variable.
	int pn, sn;                 // point number, surface number
	int rseed;                  // read in rand seed.
	double a, b, c;             // for Norm and Area Calculation.
	int p1, p2, p3;             // vertices for triangle.
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double temp;
	char path[100];     // for opening the file.
	// open the input file with the file name
	strcpy(path,fname);
	if ((fpgeo = fopen(path, "a+")) == NULL)
	{
		printf("cannot open geometry file\n");
		return 0;
	}
	// Read random seed
	fscanf(fpgeo,"%d",&rseed);
	RandSeed = (unsigned long long)rseed;
	InitRand(RandSeed);
	// read in points order(x,y,z);allocmem;readin points;Points is a 2D array.
	fscanf(fpgeo,"%d",&NumPoints);
	
	Points = (double**)malloc(NumPoints*sizeof(double*));
	for (i =  0; i < NumPoints; i++)
		Points[i] = (double*)malloc(3*sizeof(double));
	for (i = 0; i < NumPoints; i++)
	{
		fscanf(fpgeo,"%d",&pn);
		for (j = 0; j < 3; j++)
			fscanf(fpgeo, "%lf",&Points[pn][j]);
	}
	
	// read in surfaces order(p1,p2,p3,spec)
	fscanf(fpgeo,"%d",&NumSurfaces);
	dist2surf = (double*)malloc(NumSurfaces*sizeof(double));
	Surfaces = (int**)malloc(NumSurfaces*sizeof(int*));
	SurfaceMarker = (int*)malloc(NumSurfaces*sizeof(int));
	SurfaceNorms = (double**)malloc(NumSurfaces*sizeof(double*));
	SurfaceAreas = (double*)malloc(NumSurfaces*sizeof(double));
	Displacement = (double*)malloc(NumSurfaces*sizeof(double));
	Specularity = (double*)malloc(NumSurfaces*sizeof(double));
	for (i = 0; i < NumSurfaces; i++)
	{
		Surfaces[i] = (int*)malloc(3*sizeof(int));
		SurfaceNorms[i] = (double*)malloc(3*sizeof(double));
		SurfaceAreas[i] = 1.0;
		dist2surf[i] = INF;
	}
	for (i = 0; i < NumSurfaces; i++)
	{
		fscanf(fpgeo,"%d",&sn);
		for (j = 0; j < 3; j++)
		{
			fscanf(fpgeo,"%d",&Surfaces[sn][j]);
			SurfaceNorms[sn][j] = 0.0;           // set to a safe value for the norms.
		}
		fscanf(fpgeo,"%d",&SurfaceMarker[sn]);   // SurfaceMarker: mark the surface type.
		if (SurfaceMarker[sn]==0)
			Specularity[sn] = Spec;              // input outer surface specularity.
		else
			Specularity[sn] = 1.0;               // For periodic or in/out surfaces.
		// SurfaceMarker: 0: normal surfaces; 1: inflow surfaces; 2:outflow surfaces.
	}
	
	// Force define the inflow and outflow surfaces.
	fscanf(fpgeo,"%d",&NumInSurfaces);
	InSurfaces = (int*)malloc(NumInSurfaces*sizeof(int));
	for (i = 0; i < NumInSurfaces; i++)
		fscanf(fpgeo, "%d", &InSurfaces[i]);
	fscanf(fpgeo, "%d",&NumOutSurfaces);
	OutSurfaces = (int*)malloc(NumOutSurfaces*sizeof(int));
	for (i = 0; i < NumOutSurfaces; i++)
		fscanf(fpgeo, "%d", &OutSurfaces[i]);
	// read in unitcell length (in nm)
	fscanf(fpgeo, "%lf", &UnitCellLength);
	// read in the total number of particle in simulation.
	fscanf(fpgeo, "%d", &NumParticles);
	
	// First Time calculate the Surfaces Normals,Surface Area and Surface_Displacement.
	for (i = 0; i < NumSurfaces; i++)
	{
		// S = 1/2 * cross(a,b). And Norm = cross(a,b)/module.
		// get points and data
		p1 = Surfaces[i][0];
		p2 = Surfaces[i][1];
		p3 = Surfaces[i][2];
		x1 = Points[p1][0];
		y1 = Points[p1][1];
		z1 = Points[p1][2];
		x2 = Points[p2][0];
		y2 = Points[p2][1];
		z2 = Points[p2][2];
		x3 = Points[p3][0];
		y3 = Points[p3][1];
		z3 = Points[p3][2];
		// calculate the cross
		a = (y1-y2)*(z1-z3)-(z1-z2)*(y1-y3);
		b = (z1-z2)*(x1-x3)-(x1-x2)*(z1-z3);
		c = (x1-x2)*(y1-y3)-(y1-y2)*(x1-x3);
		temp = sqrt(a*a+b*b+c*c);
		
		// Nomalized Normal.
		SurfaceNorms[i][0] = a/temp;
		SurfaceNorms[i][1] = b/temp;
		SurfaceNorms[i][2] = c/temp;
		SurfaceAreas[i] = temp/2;

		// Displacement
		Displacement[i] = -(x1*SurfaceNorms[i][0]+y1*SurfaceNorms[i][1]+z1*SurfaceNorms[i][2]);
	}
	// close FILE pointer;
	fclose(fpgeo);
	// Output the data for checker;
	if ((fpchecker=fopen("./geochecker.out","a+"))==NULL)
		return -1;                 // output some error info
	fprintf(fpchecker,"-----------Output for Comparison----------\n");
	fprintf(fpchecker,"NumPoints:%d\n",NumPoints);
	for (i = 0; i<NumPoints; i++)
		fprintf(fpchecker, "%d %.6e %.6e %.6e\n",i, Points[i][0],Points[i][1],Points[i][2]);
	fprintf(fpchecker,"\n");
	fprintf(fpchecker,"NumSurfaces:%d\n",NumSurfaces);
	for (i = 0; i < NumSurfaces; i++)
		fprintf(fpchecker,"%d  %d %d %d  %.3lf  %d\n",i, Surfaces[i][0],Surfaces[i][1],Surfaces[i][2],Specularity[i],SurfaceMarker[i]);
	
	fprintf(fpchecker,"\nInflow Surfaces\n");
	fprintf(fpchecker,"NumInSurfaces:%d\n",NumInSurfaces);
	for (i = 0 ; i < NumInSurfaces; i++)
		fprintf(fpchecker, "%d\n", InSurfaces[i]);
	fprintf(fpchecker,"\nOutflow Surfaces\n");
	fprintf(fpchecker,"NumOutSurfaces:%d\n",NumOutSurfaces);
	for (i = 0; i < NumOutSurfaces; i++)
		fprintf(fpchecker, "%d\n",OutSurfaces[i]);
	fprintf(fpchecker,"\nNumPeriods:%d\n",NumPeriods);
	fprintf(fpchecker,"UnitCell Length:%.3lf\n",UnitCellLength);
	fprintf(fpchecker,"Sample Length:%.3lf\n",UnitCellLength*NumPeriods);
	fprintf(fpchecker,"Particle Number in Calculation:%d\n", NumParticles);
	for (i = 0; i < NumSurfaces; i++)
		fprintf(fpchecker,"Surface:%d, Normal=(%lf,%lf,%lf), Displacement:%lf\n",i,SurfaceNorms[i][0],SurfaceNorms[i][1],SurfaceNorms[i][2],Displacement[i]);
	fprintf(fpchecker,"-------------------END--------------------\n");
	fclose(fpchecker);
	// Build the Area Cumulative Function for inflow surfaces.
	buildCumuAreas();
	printf("Building Geometry Succeed!\n");
	// Some other operations? Output the codes to generate the 3D plot?
	/* This need only be done once. End of function*/
}


// build the cumulative function for Inflow Surfaces Area.
// Done before call transmission function, After buiding geometry
void buildCumuAreas()
{
	int i, surfnum;
	CumuAreas = (double*)malloc(NumInSurfaces*sizeof(double));
	double sum;
	sum = InSurfaces[0];
	CumuAreas[0] = sum;
	for (i = 1; i < NumInSurfaces; i++)
	{
		surfnum = InSurfaces[i];
		sum = sum + SurfaceAreas[surfnum];
		CumuAreas[i] = sum;
	}
	// Normalize
	for (i = 0; i < NumInSurfaces; i++)
		CumuAreas[i] = CumuAreas[i]/sum;
}
