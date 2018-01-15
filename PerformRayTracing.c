/* This block: perform Ray Tracing to get the transmission */
// prerequisite: Geo-builder.
// input: bulk MFP, initial position and direction.
// output: transmission rate.

#include"Functions.h"
#include"GlobalVariables.h"
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

void PerformRayTracing(double bulkMFP, particle ptc)
{
	int i,j;
	// mark the state of the particle: 0,still inside; -1,reflect; 1,transmit;
	int out = 0;                                     // mark the state of particle: in:0; transmit:1;reflect:-1
	int ZCounter=1;                                  // INIT ZCounter
	int p1,p2,p3;
	// Initialization of some of the aspects.
	double freepath;
	double cosangle;
	double S1, S2, S3;
	double dotprod;
	double c1,c2,c3,c4,s1,s2,s3,s4;                  // for diffuse scattering angle.
	double ux, uy, uz, vx, vy, vz, wx, wy, wz, u, v;       // Inside Surface or Not?
	double a, b, c, d, e;
	double theta, phi, thetaRot, phiRot;             // Angles
	double dist;                                     //
	int hitface;                                     // nearest surface.
	double position[3];                              // the new position.
	FILE *fperror;
	hitface = ptc.Surface;                      // lastly hit surface.

	// Perform Ray Tracing;
	while (!out)
	{
		// Find the next surface to collide with; one by one.
		for (i = 0; i < NumSurfaces; i++)
		{
			// calculate the angle cosine value
			a = SurfaceNorms[i][0];
			b = SurfaceNorms[i][1]; 
			c = SurfaceNorms[i][2];
			cosangle = a*ptc.tx+b*ptc.ty+c*ptc.tz;
			
			// check if the surface normal is the same direction as the particle.
			if (cosangle>0)
			{
				SurfaceNorms[i][0] = -a;
				SurfaceNorms[i][1] = -b;
				SurfaceNorms[i][2] = -c;
				Displacement[i] = -Displacement[i];
				cosangle = SurfaceNorms[i][0]*ptc.tx+SurfaceNorms[i][1]*ptc.ty+SurfaceNorms[i][2]*ptc.tz;
			}
			
			// exclude the previous surface and parallel surface
			if (i!=ptc.Surface && cosangle!=0)
			{
				// Calculate the distance to surface along the path.
				a = SurfaceNorms[i][0];
				b = SurfaceNorms[i][1];
				c = SurfaceNorms[i][2];
				// traveling distance along traveling direction.
				dist2surf[i] = -(a*ptc.x+b*ptc.y+c*ptc.z+Displacement[i])/cosangle;
				// dist2surf < 0 means it will not hit on the surface.
				if (dist2surf[i]<0)
					dist2surf[i] = INF;
				else
				{
					// Find the intercept point location.
					position[0] = ptc.x + dist2surf[i]*ptc.tx;
					position[1] = ptc.y + dist2surf[i]*ptc.ty;
					position[2] = ptc.z + dist2surf[i]*ptc.tz;
					
					// check if the intercepting point is outside the triangle.
					p1 = Surfaces[i][0];
					p2 = Surfaces[i][1];
					p3 = Surfaces[i][2];
					// U
					ux = Points[p2][0] - Points[p1][0];
					uy = Points[p2][1] - Points[p1][1];
					uz = Points[p2][2] - Points[p1][2];
					// V
					vx = Points[p3][0] - Points[p1][0];
					vy = Points[p3][1] - Points[p1][1];
					vz = Points[p3][2] - Points[p1][2];
					// W
					wx = position[0] - Points[p1][0];
					wy = position[1] - Points[p1][1];
					wz = position[2] - Points[p1][2];
					// a
					a = ux*ux + uy*uy + uz*uz;
					// b
					b = ux*vx + uy*vy + uz*vz;
					// c
					c = vx*vx + vy*vy + vz*vz;
					// d
					d = wx*ux + wy*uy + wz*uz;
					// e
					e = wx*vx + wy*vy + wz*vz;
					// u & v: coefficients
					u = (c*d-b*e)/(a*c-b*b);
					v = (a*e-b*d)/(a*c-b*b);
					//S1 = TriArea(position,Points[p1],Points[p2]);
					//S2 = TriArea(position,Points[p1],Points[p3]);
					//S3 = TriArea(position,Points[p2],Points[p3]);
					// if the intercept point is out of the triangle, it will not hit on the surface.
					//if (S1+S2+S3-SurfaceAreas[i]>1e-6)
						//dist2surf[i] = INF;
					if (u>=0 && v>=0 && u+v<=1)
						;
					else
						dist2surf[i] = INF;
				}
			}
			else
				dist2surf[i] = INF;
		}   // end of for loop: find the collided surface.

		// find the nearest surface
		dist = dist2surf[0];
		hitface = 0;
		for (i = 1; i < NumSurfaces; i++)
		{
			if (dist2surf[i]<dist)
			{
				dist = dist2surf[i];
				hitface = i;
			}
		}

		// Track Down the hitface
		//fperror = fopen("./trackparticle.txt","a+");
		//fprintf(fperror,"hitface:%d, surfacemarker:%d\n",hitface,SurfaceMarker[hitface]);
		//fclose(fperror);
		// error msg
		if (dist==INF)
		{
			fperror = fopen("./errmsg.txt","a+");
			fprintf(fperror,"Surface:%d, Particle Position(%lf,%lf,%lf),Direction(%lf,%lf,%lf)\n",ptc.Surface,ptc.x, ptc.y, ptc.z, ptc.tx, ptc.ty, ptc.tz);
			for (i = 0; i < NumSurfaces; i++)
			{
				fprintf(fperror,"Surface:#%d\n", i);
				fprintf(fperror,"Surface Normal:(%lf,%lf,%lf)\n",SurfaceNorms[i][0],SurfaceNorms[i][1],SurfaceNorms[i][2]);
				fprintf(fperror,"Surface Displacement: %lf\n", Displacement[i]);
				cosangle = SurfaceNorms[i][0]*ptc.tx+SurfaceNorms[i][1]*ptc.ty+SurfaceNorms[i][2]*ptc.tz;
				fprintf(fperror,"cosangle:%lf\n",cosangle);
				a = SurfaceNorms[i][0];
				b = SurfaceNorms[i][1];
				c = SurfaceNorms[i][2];
				dist2surf[i] = -(a*ptc.x+b*ptc.y+c*ptc.z+Displacement[i])/cosangle;
				fprintf(fperror,"dist2surf:%lf\n",dist2surf[i]);

					// intercept point
					position[0] = ptc.x + ptc.tx*dist2surf[i];
					position[1] = ptc.y + ptc.ty*dist2surf[i];
					position[2] = ptc.z + ptc.tz*dist2surf[i];
					fprintf(fperror,"intercept point:(%lf,%lf,%lf)\n",position[0],position[1],position[2]);
					p1 = Surfaces[i][0];
					p2 = Surfaces[i][1];
					p3 = Surfaces[i][2];
					S1 = TriArea(position,Points[p1],Points[p2]);
					S2 = TriArea(position,Points[p1],Points[p3]);
					S3 = TriArea(position,Points[p2],Points[p3]);
					fprintf(fperror,"S1:%lf,S2:%lf,S3:%lf, Sum:%lf, SurfaceArea:%lf\n",TriArea(position,Points[p1],Points[p2]),TriArea(position,Points[p1],Points[p3]),TriArea(position,Points[p2],Points[p3]),S1+S2+S3,SurfaceAreas[i]);
					fprintf(fperror,"position:(%lf,%lf,%lf),p1(%lf,%lf,%lf),p2(%lf,%lf,%lf),p3(%lf,%lf,%lf)\n",position[0],position[1],position[2],\
							Points[p1][0],Points[p1][1],Points[p1][2],Points[p2][0],Points[p2][1],Points[p2][2],Points[p3][0],Points[p3][1],Points[p3][2]);
			}
			fclose(fperror);
			printf("Cannot find the surface to collide with\n");
			system("pause");
		}
		// Found the surface to collide with

		// Free Path sampling
		if (bulkMFP==INF)
			freepath = INF;
		else
			freepath = -bulkMFP*log(RandR());
		
		// Calculate the new position;
		ptc.x = ptc.x + ptc.tx*dist;
		ptc.y = ptc.y + ptc.ty*dist;
		ptc.z = ptc.z + ptc.tz*dist;
		
		// intrinsic scattering
		if (freepath < dist)
		{
			dist = freepath-dist;
			// updating the position.
			ptc.x = ptc.x + ptc.tx*dist;
			ptc.y = ptc.y + ptc.ty*dist;
			ptc.z = ptc.z + ptc.tz*dist;
			
			// randomly set the new traveling direction
			theta = acos(1-2*RandR());              // Math functions like asin and sqrt need to be defined.
			phi = 2*PI*RandR();                       // PI need to be defined.
			ptc.tx = sin(theta)*cos(phi);
			ptc.ty = sin(theta)*sin(phi);
			ptc.tz = cos(theta);                      // update traveling direction.
			ptc.Surface = -1;                         // set the last hit surface as -1;
		}
		// collide with surfaces.
		else if (SurfaceMarker[hitface]==1)            // SurfaceMarker=1: In-flow Surfaces.
			{
				if (ZCounter==1)                      // exit from InSurfaces
					out = -1;
				else
				{
					ZCounter = ZCounter-1;            // go back the previous unit cell
					// traveling direction doesn't change;
					// position @ Z move to the length.
					ptc.z = UnitCellLength;       // move z to the end of the simulation box.
					// Here assume the traveling direction is along +Z direction.
					// Find the hit surfaces from the out-flow surfaces: use j
					for (j = 0; j < NumOutSurfaces; j++)
					{
						if (IsInsideTriangle2D(ptc.x, ptc.y, OutSurfaces[j]))
						{
							//printf("Last Hit:%d, Fake Hit:%d\n",ptc.Surface, hitface);
							ptc.Surface = OutSurfaces[j];
							//printf("Found the OutSurfaces:%d!\n",ptc.Surface);
							break;
						}
					}
					/*if (!Found)
					{
						//printf("Didn't Found the Corresponding Surface In OutSurfaces!\n");
						//system("pause");
					}*/
				}
			}
		// type II
		else if (SurfaceMarker[hitface]==2)            // SurfaceMarker=2: Out-flow Surfaces.
			{
				if (ZCounter==NumPeriods)             // exit from the OutSurfaces
					out = 1;
				else
				{
					ZCounter = ZCounter + 1;          // Move to the next unit cell
					// traveling direction doesn't change
					// position @ Z move to the 0;
					ptc.z = 0;
					// here assume traveling direction is along +Z direction.
					// find the corresponding surface from In-flow Surfaces. use j;
					for (j = 0; j < NumInSurfaces; j++)
					{
						if (IsInsideTriangle2D(ptc.x,ptc.y,InSurfaces[j]))
						{
							//printf("Last Hit:%d, Fake Hit:%d\n",ptc.Surface, hitface);
							ptc.Surface = InSurfaces[j];
							//printf("Found the OutSurfaces:%d!\n",ptc.Surface);
							break;
						}
					}
					/*if (!Found)
					{
						//printf("Didn't find the corresponding surfaces in InSurfaces\n");
						///system("pause");
					}*/
				}
			}
			// 	Regular surfaces.
		else
		{
			ptc.Surface = hitface;
			// specularity; BS rate included or not?
			// specular
			if (RandR()<=Specularity[hitface])
			{
				// change direction.
				dotprod = -2*(ptc.tx*SurfaceNorms[hitface][0]+ptc.ty*SurfaceNorms[hitface][1]+ptc.tz*SurfaceNorms[hitface][2]);
				ptc.tx += dotprod*SurfaceNorms[hitface][0];
				ptc.ty += dotprod*SurfaceNorms[hitface][1];
				ptc.tz += dotprod*SurfaceNorms[hitface][2];
				// Nothing else.
			}
			// diffuse
			else
			{
				// Firstly generate the direction and the multiply rotation matrix.
				// original angles.
				theta = asin(sqrt(RandR()));                  //
				phi = 2*PI*RandR();                           // PI is constant.
				// X,Y,Z partial directions in original Coordinate.
				// Rotation Angles
				phiRot = atan2(SurfaceNorms[hitface][1],SurfaceNorms[hitface][0]);    // thetaRot = atan2(ny,nx); atan2 to be defined.
				thetaRot = acos(SurfaceNorms[hitface][2]);                                // phiRot = asin(nz);
				// the cosine and sine value for these four angles
				s1 = sin(theta);
				c1 = cos(theta);
				s2 = sin(phi);
				c2 = cos(phi);
				s3 = sin(thetaRot);
				c3 = cos(thetaRot);
				s4 = sin(phiRot);
				c4 = cos(phiRot);
				// Directly have the direction after reflection.
				ptc.tx = (s1*c2*c3*c4-s1*s2*s4+c1*s3*c4);
				ptc.ty = (s1*c2*c3*s4+s1*s2*c4+c1*s3*s4);
				ptc.tz = (-s1*c2*s3+c1*c3);
				// hit the surface, not need to update the position.
			}
		}
	}

	// Outside the while loop
	if (out==1)
	{
		//printf("Transmit\n");
		NumTrans = NumTrans + 1;
	}
	else
	{
		//printf("Reflect\n");
		NumReflect = NumReflect + 1;
	}
	// There seem to be nothing we can do.
}

// IsInsideTriangle2D
// Input: position of the point, the surface number.
// Output: true:1; false:0
int IsInsideTriangle2D(double px, double py, int SurfaceNum)
{
	int p1, p2, p3;               // three vertices of the triangle.
	int n = SurfaceNum;           // surface number.
	double a1,a2,b1,b2,c1,c2;     // for telling the inside or not
	double t1, t2, t3;
	p1 = Surfaces[n][0];
	p2 = Surfaces[n][1];
	p3 = Surfaces[n][2];
	a1 = px - Points[p1][0];
	a2 = py - Points[p1][1];
	b1 = px - Points[p2][0];
	b2 = py - Points[p2][1];
	c1 = px - Points[p3][0];
	c2 = py - Points[p3][1];
	t1 = (a1*b2-a2*b1);
	t2 = (b1*c2-b2*c1);
	t3 = (c1*a2-c2*a1);
	if (t1*t2>0 && t1*t3>0)
		return 1;
	else
		return 0;
}
