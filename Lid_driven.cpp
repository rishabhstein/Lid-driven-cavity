#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define grid 50

int main (void)
{
	double u[grid][grid+1], un[grid][grid+1], uc[grid][grid];
	double v[grid+1][grid], vn[grid+1][grid], vc[grid][grid];
	double p[grid+1][grid+1], pn[grid+1][grid+1], pc[grid][grid];
	double m[grid+1][grid+1];
	int i, j, step;
	double dx, dy, dt, tau, delta, error, Re;
	step =1;
	dx = 1.0/(grid-1);
	dy = 1.0/(grid-1);
	dt = 0.001;
	delta = 4.5;
	error = 1.0;
	Re = 10.0;
	
	// Initializing u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				u[i][j] = 0.0;
				u[i][grid] = 1.0;
				u[i][grid-1] = 1.0;
			}
		}
		
	// Initializing v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				v[i][j] = 0.0;
			}
		}
		
	// Initializing p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				p[i][j] = 1.0;
			}
		}
	
	while (error > 0.000001)
	{
		// Solve u-momentum equation
		for (i=1; i<=(grid-2); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				un[i][j] = u[i][j] - dt*(  (u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/2.0/dx 
							+0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy  )
								- dt/dx*(p[i+1][j]-p[i][j]) 
									+ dt*1.0/Re*( (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy );
			}
		}
		
		// Boundary conditions
		for (j=1; j<=(grid-1); j++)
		{
			un[0][j] = 0.0;
			un[grid-1][j] = 0.0;
		}
		
		for (i=0; i<=(grid-1); i++)
		{
			un[i][0] = -un[i][1];
			un[i][grid] = 2 - un[i][grid-1];
		}
		
		
		// Solves v-momentum
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-2); j++)
			{
				vn[i][j] = v[i][j] - dt* ( 0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j]) )/dx 
							+(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/2.0/dy ) 
								- dt/dy*(p[i][j+1]-p[i][j]) 
									+ dt*1.0/Re*( (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy );
			}
		}
		
		// Boundary conditions
		for (j=1; j<=(grid-2); j++)
		{
			vn[0][j] = -vn[1][j];
			vn[grid][j] = -vn[grid-1][j];
		}		

		for (i=0; i<=(grid); i++)
		{
			vn[i][0] = 0.0;
			vn[i][grid-1] = 0.0;
		}		
	
		// Solves continuity equation
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				pn[i][j] = p[i][j]-dt*delta*(  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] ) /dy  );
			}
		}
		
		
		// Boundary conditions
		for (i=1; i<=(grid-1); i++)
		{
			pn[i][0] = pn[i][1];
			pn[i][grid] = pn[i][grid-1];
		}
		
		for (j=0; j<=(grid); j++)
		{
			pn[0][j] = pn[1][j];
			pn[grid][j] = pn[grid-1][j];
		}		
		
		// Displaying error
		error = 0.0;
		
		for (i=1; i<=(grid-1); i++)
		{
			for (j=1; j<=(grid-1); j++)
			{
				m[i][j] = (  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] )/dy  );
				error = error + fabs(m[i][j]);
			}
		}
		
		if (step%10000 ==1)
		{
	    printf("Error is %5.8lf for the step %d\n", error, step);
		}
		
		
		// Iterating u
		for (i=0; i<=(grid-1); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				u[i][j] = un[i][j];
			}
		}
		
		// Iterating v
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid-1); j++)
			{
				v[i][j] = vn[i][j];
			}
		}
		
		// Iterating p
		for (i=0; i<=(grid); i++)
		{
			for (j=0; j<=(grid); j++)
			{
				p[i][j] = pn[i][j];
			}
		}

		step = step + 1;
	
	}
	
	for (i=0; i<=(grid-1); i++)
	{
		for (j=0; j<=(grid-1); j++)
		{	
			uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
			vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
			pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		}
	}
	
	
	
	// OUTPUT DATA
	FILE *fout2, *fout3 ,*fout4;
	fout2 = fopen("UVP.plt","w+t");
	fout3 = fopen("Central_U.plt","w+t");
	fout4 = fopen("Central_V.plt","w+t");

	if ( fout2 == NULL )
	{
    printf("\nERROR when opening file\n");
    fclose( fout2 );
	}

  else
	{
	fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	fprintf( fout2, "ZONE  F=POINT\n");
	fprintf( fout2, "I=%d, J=%d\n", grid, grid );

	for ( j = 0 ; j < (grid) ; j++ )
	{
    for ( i = 0 ; i < (grid) ; i++ )
    {
		double xpos, ypos;
		xpos = i*dx;
		ypos = j*dy;

		fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uc[i][j], vc[i][j], pc[i][j] );
    }
	}
	}

	fclose( fout2 );
	
	// CENTRAL --U
  fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
  fprintf(fout3, "ZONE F=POINT\n");
  fprintf(fout3, "I=%d\n", grid );

  for ( j = 0 ; j < grid ; j++ )
  {
	double ypos;
    ypos = (double) j*dy;

    fprintf( fout3, "%5.8lf\t%5.8lf\n", (uc[grid/2][j] + uc[(grid/2)+1][j])/(2.), ypos );
  }
  //CENTRAL --V
  fprintf(fout4, "VARIABLES=\"V\",\"X\"\n");
  fprintf(fout4, "ZONE F=POINT\n");
  fprintf(fout4, "J=%d\n", grid );

  for ( i = 0 ; i < grid ; i++ )
  {
	double xpos;
    xpos = (double) i*dy;

    fprintf( fout4, "%5.8lf\t%5.8lf\n", (vc[i][grid/2] + vc[i][(grid/2)+1])/(2.), xpos );
  }
}
