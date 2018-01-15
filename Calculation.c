/* This block of code read in bulk thermal conductivity MFP spectra and
launch the Ray Tracing to calculate the effective MFP for nanostructure
, and complish the integration procedure to get the thermal conductivity*/

// input file: bulktcspectra.in
// output file: Effective BS rate plot.
#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include"GlobalVariables.h"
#include"Functions.h"
#include<string.h>
#include<math.h>

#define MAXFILENAME 100

int CalTherm(char *fname)
{
	FILE *fpspectra;                      // read in the spectra
	FILE *fpchecker;                      // output the readin info to check code.
	FILE *fptc;                           // output the results for thermal conductivity
	FILE *fpresult;                        // output the transmission results.
	int i;                                // loop variable
	double bulkMFP, effectiveMFP, transmission;
	double delta;
	double K;
	double SampleLength;
	time_t start, finish;
	char path[MAXFILENAME];
	SampleLength = NumPeriods*UnitCellLength;

	// get clock
	start = clock();
	strcpy(path,fname);
	// read in the TC-MFP spectra
	if ((fpspectra=fopen(path,"r"))==NULL)
		return -1;
	/* TC spectra form: 
	MFP K_Lambda dMFP
	*/
	fscanf(fpspectra,"%d", &NumBins);
	MFPs = (double*)malloc(NumBins*sizeof(double));
	Klambda = (double*)malloc(NumBins*sizeof(double));
	dlambda = (double*)malloc(NumBins*sizeof(double));
	for (i = 0; i < NumBins; i++)
	{
		fscanf(fpspectra, "%lf",&MFPs[i]);
		fscanf(fpspectra, "%lf",&Klambda[i]);
		fscanf(fpspectra, "%lf",&dlambda[i]);
	}
	fclose(fpspectra);
	
	// output the readin info to check
	if ((fpchecker=fopen("./bulktcspectra.out","a+"))==NULL)
		return -1;
	fprintf(fpchecker,"---------------TC_spectra-----------------\n");
	fprintf(fpchecker,"Number of MFP bins: %d\n",NumBins);
	for (i = 0; i < NumBins; i++)
	{
		fprintf(fpchecker, "%.4e  ",MFPs[i]);
		fprintf(fpchecker, "%.4e  ", Klambda[i]);
		fprintf(fpchecker, "%.4e\n", dlambda[i]);
	}
	fprintf(fpchecker,"\n+++++++++++++++++End_of_File++++++++++++++++\n");
	fclose(fpchecker);
	
	// Call transmission calculation function & conduct numerical integration.
	effectiveTC = 0;
	if ((fpresult=fopen("./results.dat","a+"))==NULL)
		return 0;          // and other error msg.
	fprintf(fpresult,"Variables = \"Bulk_MFP\",\"Nano_MFP\",\"Transmission\"\n\n");
	fprintf(fpresult,"Zone T = \"RayTracing\", I = %d, DataPacking = Point\n", NumBins);
	for (i = 0; i < NumBins; i++)
	{
		K = Klambda[i];
		bulkMFP = MFPs[i];
		delta = dlambda[i];
		transmission = CalTransmission(bulkMFP);
		effectiveMFP = 3.0/4.0*transmission*SampleLength;
		effectiveTC = effectiveTC + K*delta*effectiveMFP/bulkMFP;
		fprintf(fpresult,"%lf  %lf   %lf\n",bulkMFP, effectiveMFP, transmission);
	}
	fclose(fpresult);
	// output the final result
	if ((fptc=fopen("./ThermalConductivity.out","a+"))==NULL)
		return 0;   // and other error msg.
	finish = clock();
	fprintf(fptc,"Effective Thermal Conductivity: %.6e[W/mK]\n",effectiveTC);
	fprintf(fptc,"Total Time Elapsed: %.3lfs\n",(double)(finish-start)/CLOCKS_PER_SEC);
	fclose(fptc);
	return 0;
}

void CalMFPs(char *fname)
{
	FILE *fpmfps;                                // Read in MFPs.
	FILE *fpchecker;                             // Output readin info to check.
	FILE *fpresult;                                 // Output results.
	int i;
	double effectiveMFP, bulkMFP,transmission;
	double SampleLength;
	char path[MAXFILENAME];
	SampleLength = UnitCellLength*NumPeriods;

	time_t start, finish;
	start = clock();
	strcpy(path,fname);
	// read in the infos
	if ((fpmfps=fopen(path,"r"))==NULL)
	{
		printf("Error in opening MFP files\n");
	}
	fscanf(fpmfps,"%d",&NumBins);
	MFPs = (double*)malloc(NumBins*sizeof(double));
	for (i = 0; i < NumBins; i++)
		fscanf(fpmfps, "%lf", &MFPs[i]);
	fclose(fpmfps);
	
	// output info to check
	if ((fpchecker=fopen("./MFPchecker.out","a+"))==NULL)
	{
		printf("Error in creating MFPchecker.out\n");
	}
	fprintf(fpchecker,"\n---------------Input_MFP-----------------\n");
	fprintf(fpchecker,"NumBins:%d\n",NumBins);
	for (i = 0; i < NumBins; i++)
		fprintf(fpchecker,"%lf\n",MFPs[i]);
	fprintf(fpchecker,"-------------------EOF---------------------\n");
	fclose(fpchecker);
	// Calculate effective MFPs one by one.
	if ((fpresult=fopen("./Results.dat","a+"))==NULL)
	{
		printf("Error in creating Results.out\n");
	}
	/*
	output form to Tecplot
	bulkMFP  Transmission  EffectiveMFP
	*/
	fprintf(fpresult, "Variables = \"Bulk_MFP\",\"Nano_MFP\",\"Transmission\"\n\n");
	fprintf(fpresult, "Zone T = \"RayTracing\", I = %d, DataPacking = Point\n", NumBins);
	for (i = 0; i < NumBins; i++)
	{
		bulkMFP = MFPs[i];
		transmission = CalTransmission(bulkMFP);
		effectiveMFP = 3.0/4.0*transmission*SampleLength;
		fprintf(fpresult,"%.3e  %.3e   %.3lf\n",bulkMFP, effectiveMFP, transmission);
	}
	finish = clock();
	fprintf(fpresult,"time elapsed:%lf s\n", (double)(finish-start)/CLOCKS_PER_SEC);
	fclose(fpresult);
}

// Calculate only the bdy MFP for length from 1 to 20 Period
void CalBdyMFPs(char *pathgeo)
{
	FILE *fpresults, *fpmfps;
	int i, j, NumBins;
	double *MFPs;
	double transmission,effectiveMFP,bulkMFP,SampleLength;
	if ((fpmfps = fopen("./mfps.in","r"))==NULL)
	{
		printf("Cannot open ./mfps.in\n");
		exit(0);
	}
	fscanf(fpmfps,"%d", &NumBins);
	MFPs = (double*)malloc(NumBins*sizeof(double));
	for (i = 0; i < NumBins; i++)
		fscanf(fpmfps,"%lf",&MFPs[i]);
	fclose(fpmfps);
	// checker mfps
	fpmfps = fopen("./CheckMFPs.txt","a+");
	fprintf(fpmfps,"NumBins:%d\n",NumBins);
	for (i = 0; i < NumBins; i++)
		fprintf(fpmfps,"%.3lf\n",MFPs[i]);
	fclose(fpmfps);
	Spec = 0;
	fpresults = fopen("./OutMFPs.txt","a+");
	fprintf(fpresults,"Specularity=0\n| NumPeriods | bulk_MFP | MFP_bdy | Transmission |\n");
	fclose(fpresults);
	for (i = 1; i <= 30; i++)
	{
		NumPeriods = i;
		GeoBuilder(pathgeo);
		for (j = 0; j < NumBins; j++)
		{
			fpresults = fopen("./OutMFPs.txt","a+");
			bulkMFP = MFPs[j];
			SampleLength = UnitCellLength*NumPeriods;
			transmission = CalTransmission(bulkMFP);
			effectiveMFP = 3.0/4.0*transmission*SampleLength;
			fprintf(fpresults,"%d  %.3lf  %e   %.3f\n", i, bulkMFP, effectiveMFP, transmission);
			fclose(fpresults);
		}
	}
	/*
	Spec = 1;
	fpresults = fopen("./OutMFPs.txt","a+");
	fprintf(fpresults,"Specularity=1\n| NumPeriods | bulk_MFP | MFP_bdy | Transmission |\n");
	fclose(fpresults);
	for (i = 1; i <= 50; i++)
	{
		NumPeriods = i;
		GeoBuilder(pathgeo);
		for (j = 0; j < NumBins; j++)
		{
			fpresults = fopen("./OutMFPs.txt","a+");
			bulkMFP = INF;
			SampleLength = UnitCellLength*NumPeriods;
			transmission = CalTransmission(bulkMFP);
			effectiveMFP = 3.0/4.0*transmission*SampleLength;
			fprintf(fpresults,"%d  %.3lf  %e   %.3f\n", i, bulkMFP, effectiveMFP, transmission);
			fclose(fpresults);
		}
	}*/
}
// Calculate triangle areas
double TriArea(double p1[3], double p2[3], double p3[3])
{
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double a, b, c, area;
	x1 = p1[0]; y1 = p1[1]; z1 = p1[2];
	x2 = p2[0]; y2 = p2[1]; z2 = p2[2];
	x3 = p3[0]; y3 = p3[1]; z3 = p3[2];
	a = (y1-y2)*(z1-z3)-(z1-z2)*(y1-y3);
	b = (z1-z2)*(x1-x3)-(x1-x2)*(z1-z3);
	c = (x1-x2)*(y1-y3)-(y1-y2)*(x1-x3);
	area = sqrt(a*a+b*b+c*c);
	return (area/2);
}
