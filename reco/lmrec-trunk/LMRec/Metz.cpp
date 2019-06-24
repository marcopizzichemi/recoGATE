/*
  Name: Metz.cpp
  Copyright: 
  Author: Liji Cao
  Date: 03/11/11 11:06
  Description: These functions perform Metz filtering with given sigma of Gaussian and Metz power. 
               Images are automatically converted to FFTW COMPLEX type, filtered, and converted back to the original type.
*/

#include<math.h>
#include <fftw3.h>
#include<stdlib.h>
#include<stdio.h>

#define pi M_PI

void Metz2D(float * Data, int Xdim, int Ydim, float StdDev, float Power)
{
   float kx = 2*pi / (float)Xdim;
   float ky = 2*pi / (float)Ydim;

   int x, y, x1, y1, halfx, halfy;
   double numerator, denominator, H, Q;
   float total = 0;
   for (int i=0; i<Xdim; i++)
        for (int j=0; j<Ydim; j++)
            total += Data[i+j*Xdim];

   fftw_complex *fftwdata = (fftw_complex*)fftw_malloc( Xdim*Ydim*sizeof(fftw_complex) );
   for (int i=0; i<Xdim; i++)
        for (int j=0; j<Ydim; j++)
        {
            fftwdata[i+j*Xdim][0] = Data[i+j*Xdim];
            fftwdata[i+j*Xdim][1] = 0;
        }
    fftw_complex *fftwout = (fftw_complex*)fftw_malloc( Xdim*Ydim*sizeof(fftw_complex) );

    fftw_plan plan;
    plan = fftw_plan_dft_2d(Xdim, Ydim, fftwdata, fftwout, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

   halfx = Xdim / 2;
   halfy = Ydim / 2;

   /* Loop over Y axis */
   for (y = 0; y < Ydim; y++)
   {
      if (y < halfy)
	 y1 = y;
      else
	 y1 = Ydim - y;

      /* Loop over X axis */
      for (x = 0; x < Xdim; x++)
      {
	 if (x < halfx)
	    x1 = x;
	 else
	    x1 = Xdim - x;

	 /* Calculate value of unit height Gaussian H */
	 numerator = kx * kx * x1 * x1 + ky * ky * y1 * y1;
	 H = exp(- numerator * StdDev * StdDev / 2.0);

	 /* Do pointwise multiplication */
	 if (H == 0)
	    Q = 0;
	 else
	    Q = (1 - pow((double) (1 - H * H), (double) (Power + 1))) / H;

	 fftwout[x+y*Xdim][0] *= Q;
	 fftwout[x+y*Xdim][1] *= Q;
      }
   }

   plan = fftw_plan_dft_2d(Xdim, Ydim, fftwout, fftwdata, FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(plan);

   for (int i=0; i<Xdim; i++)
        for (int j=0; j<Ydim; j++)
            Data[i+j*Xdim] = sqrt(fftwdata[i+j*Xdim][0]*fftwdata[i+j*Xdim][0]+fftwdata[i+j*Xdim][1]*fftwdata[i+j*Xdim][1]);

    float newTotal = 0;
    for (int i=0; i<Xdim; i++)
        for (int j=0; j<Ydim; j++)
            newTotal += Data[i+j*Xdim];

    for (int i=0; i<Xdim; i++)
        for (int j=0; j<Ydim; j++)
            Data[i+j*Xdim] *= total/newTotal;
}

void Metz3D(float * Data, int Xdim, int Ydim, int Zdim, float StdDev, float Power)
{
   float kx = 2*pi / (float)Xdim;
   float ky = 2*pi / (float)Ydim;
   float kz = 2*pi / (float)Zdim;

   int x, y, z, x1, y1, z1, halfx, halfy, halfz;
   double numerator, denominator, H, Q;
   float total = 0;
   for (int i=0; i<Xdim*Ydim*Zdim; i++)
       total += Data[i];

   fftw_complex *fftwdata = (fftw_complex*)fftw_malloc( Xdim*Ydim*Zdim*sizeof(fftw_complex) );

   for (int i = 0; i < Xdim; i++)
      for (int j = 0; j < Ydim; j++)
         for (int k = 0; k < Zdim; k++) {
            fftwdata[k+j*Zdim+i*Zdim*Ydim][0] = Data[i+j*Xdim+k*Xdim*Ydim];
            fftwdata[k+j*Zdim+i*Zdim*Ydim][1] = 0;
         }

    fftw_complex *fftwout = (fftw_complex*)fftw_malloc( Xdim*Ydim*Zdim*sizeof(fftw_complex) );

    fftw_plan plan;
    plan = fftw_plan_dft_3d(Xdim, Ydim, Zdim, fftwdata, fftwout, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

   halfx = Xdim / 2;
   halfy = Ydim / 2;
   halfz = Zdim / 2;

   for (z = 0; z < Zdim; z++)
   {
       if (z < halfz)
          z1 = z;
       else
           z1 = Zdim - z;

       for (y = 0; y < Ydim; y++)
       {
           if (y < halfy)
	          y1 = y;
           else
	           y1 = Ydim - y;

           for (x = 0; x < Xdim; x++)
           {
	           if (x < halfx)
	              x1 = x;
               else
	               x1 = Xdim - x;

               /* Calculate value of unit height Gaussian H */
	           numerator = kx * kx * x1 * x1 + ky * ky * y1 * y1 + kz * kz * z1 * z1;
	           H = exp(- numerator * StdDev * StdDev / 2.0);

	           /* Do pointwise multiplication */
	           if (H == 0)
	              Q = 0;
               else
	               Q = (1 - pow((double) (1 - H * H), (double) (Power + 1))) / H;

               fftwout[z+y*Zdim+x*Zdim*Ydim][0] *= Q;
	       fftwout[z+y*Zdim+x*Zdim*Ydim][1] *= Q;
           }
       }
   }

   plan = fftw_plan_dft_3d(Xdim, Ydim, Zdim, fftwout, fftwdata, FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(plan);

   for (int i=0; i<Xdim; i++)
        for (int j=0; j<Ydim; j++)
            for (int k=0; k<Zdim; k++) {
		Data[i+j*Xdim+k*Xdim*Ydim] = sqrt(
                   fftwdata[k+j*Zdim+i*Zdim*Ydim][0] * fftwdata[k+j*Zdim+i*Zdim*Ydim][0] +
                   fftwdata[k+j*Zdim+i*Zdim*Ydim][1] * fftwdata[k+j*Zdim+i*Zdim*Ydim][1]
                  );
            }

    float newTotal = 0;
    for (int i=0; i<Xdim*Ydim*Zdim; i++)
                newTotal += Data[i];
    for (int i=0; i<Xdim*Ydim*Zdim; i++)
                Data[i] *= total/newTotal;
}
