/* Free Memory */
// End of the program: cleanup the memory space
// Input:
// Output None

#include<stdio.h>
#include<stdlib.h>
#include"Functions.h"
#include"GlobalVariables.h"

void FreeMem(int Indicator)
{
	//inidicator = 0: doesn't calculate TC; indicator = 1: Calculate TC directly.
	int i;
	for (i = 0; i < NumPoints; i++)
		free(Points[i]);
	free(Points);
	for (i = 0; i < NumSurfaces; i++)
	{
		free(SurfaceNorms[i]);
		free(Surfaces[i]);
	}
	free(SurfaceNorms);
	free(Surfaces);
	free(SurfaceAreas);
	free(Specularity);
	free(InSurfaces);
	free(OutSurfaces);
	free(Displacement);
	free(SurfaceMarker);
	free(CumuAreas);
	free(MFPs);
	free(dist2surf);
	// extra work to be done.
	if (Indicator == 1)
	{
		free(Klambda);
		free(dlambda);
	}
}
