#include<stdio.h>
#include<stdlib.h>
#include"GlobalVariables.h"
#include"Functions.h"

int main(void)//int argc, char* argv[])
{
	// depends on the parameters input to decide what to do.
	// input: ./RayTracing indicator(Cal TC or not), geofile,  MFPs(or spectrum).
	char *pathgeo, *pathspectra;
	int indicator;
	indicator = 0;//atoi(argv[1]);
	// Does not calculate thermal conductivity
	if (indicator == 0)
	{
		pathgeo = "geo.in";
		pathspectra = "mfps.in";
		Spec = 0;
		GeoBuilder(pathgeo);
		printf("Succeeded Building Geometry\n");
		// Conduct Ray Tracing only once:
		CalMFPs(pathspectra);
		FreeMem(0);
	}
	// Calculate the thermal conductivity directly
	else if (indicator== 1)
	{
		// Firstly build the geometry
		pathgeo = "./geo.in";// argv[2];
		pathspectra = "Si.in";// = argv[3];
		//NumPeriods = atoi(argv[4]);
    	//Spec = atoi(argv[5]);
		GeoBuilder(pathgeo);
		CalTherm(pathspectra);
		FreeMem(1);
	}
	// error msg: wrong inputs
	else if (indicator == 2)
	{
		pathgeo = "geo.in";//argv[2];
		CalBdyMFPs(pathgeo);
		FreeMem(0);
	}
	else
	{
		printf("Error in input parameters:\n");
		printf("Correct input form: ./function (0:Only MFP;1:TC) in.geometry in.mfps/in.spectrum\n");
		return 0;
	}
	printf("Calculation Done!\n");
	return 1;
}
