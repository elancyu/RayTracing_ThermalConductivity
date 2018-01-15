/* This file initialize transmission calculations. */
#include"Functions.h"
#include<stdio.h>
#include<math.h>
#include"GlobalVariables.h"

double CalTransmission(double bulkMFP)
{
	int i, j, insurf, surfnum;                  
	int p1, p2, p3;                             // three vertices of triangle surface.
	double ca, cb, cc;                          // coefficients.
	double R;                                   // Random number for choosing the in-flow surface.
	double R1, R2;                              // Random number for choosing position.
	double theta, phi;                          // for traveling direction.
	double transmission;                        // return value; transmission rate.
	double costheta, sintheta, cosphi, sinphi;  // sin(theta), cos(theta), sin(phi), cos(phi)
	int fracnum;
	particle ptc;
	// Initialization: set the NumTrans and NumReflect as 0;
	NumTrans = 0;
	NumReflect = 0;
	fracnum = NumParticles*5/100;              // 5 percent

	// Calculate transmission by performing Ray Tracing
	insurf = 0;                     // default as the first surface; then search to update it.
	for (i = 0; i < NumParticles; i++)
	{
		if ((i+1)%fracnum == 0)
			printf("%d %% is compeleted\n",i/fracnum*5);
		// firstly choose which surface it is on
		R = RandR();
		for (j = NumInSurfaces; j > 0; j--)
			if (R>=CumuAreas[j-1])
			{
				insurf = j;
				break;                           // stop searching.
			}
		// The surface is InSurfaces[insurf];
		R1 = RandR();
		R2 = RandR();
		surfnum = InSurfaces[insurf];
		ptc.Surface = surfnum;
		p1 = Surfaces[surfnum][0];
		p2 = Surfaces[surfnum][1];
		p3 = Surfaces[surfnum][2];
		// coefficients
		ca = 1 - sqrt(R1);
		cb = sqrt(R1)*(1-R2);
		cc = R2*sqrt(R1);
		// Initialize the Positions.
		ptc.x = ca*Points[p1][0] + cb*Points[p2][0] + cc*Points[p3][0];
		ptc.y = ca*Points[p1][1] + cb*Points[p2][1] + cc*Points[p3][1];
		ptc.z = ca*Points[p1][2] + cb*Points[p2][2] + cc*Points[p3][2];
		// Initialize traveling direction:
		theta = asin(sqrt(RandR()));
		phi = 2*PI*RandR();
		sintheta = sin(theta);
		costheta = cos(theta);
		sinphi = sin(phi);
		cosphi = cos(phi);
		ptc.tx = sintheta*cosphi;
		ptc.ty = sintheta*sinphi;
		ptc.tz = costheta;
		PerformRayTracing(bulkMFP, ptc);
		//printf("%d is out\n",i);
	}
	printf("bulkMFP:%lf, NumTrans:%d, NumReflect:%d, NumParticles:%d\n", bulkMFP, NumTrans, NumReflect, NumParticles);
	transmission = NumTrans*1.0/NumParticles;
	return (transmission);
}
