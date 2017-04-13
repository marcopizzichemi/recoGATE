#include "lmrec_common.hpp"
#include "pem_geometry.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <string.h>

void * rec_kernel(void* parameter)
{
	struct rec_kernel_parameter *data = (struct rec_kernel_parameter*)parameter;

	int nRays = data->nRays;
	boost::mt19937 *generator = data->generator;
	boost::uniform_real<> range(-1, 1);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > next(*generator, range);


	float centerX = PEM_xLength/2;
	float centerY = PEM_yLength/2;

	srand( (unsigned)time( NULL ) ); 

	float randomDriftX1 = 0, randomDriftY1 = 0, randomDriftZ1 = 0;
	float randomDriftX2 = 0, randomDriftY2 = 0, randomDriftZ2 = 0;

	for (int i=data->beginNum; i<data->totalNum + data->beginNum; i++)	
	{

		DKFZFormat &event = data->listData[i];

		float cosTheta = cosf(event.angle);
		float sinTheta = sinf(event.angle);


		if(nRays <= 0) {
			float dx = 0;
			float dy = 0;
		
			float x1 = (event.x1+dx); 
			float y1 = (event.y1+dy) * cosTheta + event.z1 * sinTheta;
			float z1 = - (event.y1+dy) * sinTheta + event.z1 * cosTheta;
			
			float x2 = (event.x2+dx);
			float y2 = (event.y2+dy) * cosTheta + event.z2 * sinTheta;
			float z2 = - (event.y2+dy) * sinTheta + event.z2 * cosTheta;
			
			siddon_rec1( data->ip,
				     x1, y1, z1,
				     x2, y2, z2,
				     data->recImage,data->diffImage,
				     event.weight,
					 data->randImage,
					 data->rand_cut
				);
		}
		else {
			for(int ray = 0; ray < nRays; ray++) {
				float dx = next();
				float dy = next();

				float x1 = (event.x1+dx); 
				float y1 = (event.y1+dy) * cosTheta + event.z1 * sinTheta;
				float z1 = - (event.y1+dy) * sinTheta + event.z1 * cosTheta;
				
				float x2 = (event.x2+dx);
				float y2 = (event.y2+dy) * cosTheta + event.z2 * sinTheta;
				float z2 = - (event.y2+dy) * sinTheta + event.z2 * cosTheta;
				
				siddon_rec1( data->ip,
					     x1, y1, z1,
					     x2, y2, z2,
					     data->recImage,data->diffImage,
					     event.weight/nRays,
						 data->randImage,
						 data->rand_cut
					);
			}
			
		}
		
	}
	pthread_exit(NULL);
}



void siddon_rec1(struct image_params *ip, 
		 float x1, float y1, float z1, float x2, float y2, float z2, 
		 float *image, float *diffImage, float norm, float *randImage, float rand_cut)
{
	float PIXELLENGTH = ip->pixelLength;
	int XDim = ip->xDim;
	int YDim = ip->yDim;
	int ZDim = ip->zDim;

	x1 /= float(PIXELLENGTH);
	y1 /= float(PIXELLENGTH);
	z1 /= float(PIXELLENGTH);
	x2 /= float(PIXELLENGTH);
	y2 /= float(PIXELLENGTH);
	z2 /= float(PIXELLENGTH);
	register int x, y, z = 0, bv;
	int   bilstartx, bilstarty, bilstartz, bilstopx, bilstopy, bilstopz,
		BV, hfsx = XDim / 2, hfsy = YDim / 2, hfsz = ZDim / 2;
	char  xflag, yflag, zflag;
	float ax[XDim*2], ay[YDim*2], az[ZDim*2], bx, by, bz, blpx, nsux, 
		nsuy, nsuz, psux, psuy, psuz, dmin, dminx, dminy, dminz, dmax, 
		dmaxx, dmaxy, dmaxz, tx, ty, tz, batte, gammatogammaforw, amid, 
		blivcm, blivpx, bvattecoeff, battefac, bvattew_, fact, t1, t2;
	int          *x_ = NULL, *y_ = NULL, *z_ = NULL;
	float        *alpha = NULL, *bvattew = NULL;

	gammatogammaforw = 0.0;
	float randomforw = 0.0;
	float randnorm;
 
	if (fabs(bx = x1 - x2) < 1.e-6) bx = 1.e-6;
	if (fabs(by = y1 - y2) < 1.e-6) by = 1.e-6;
	if (fabs(bz = z1 - z2) < 1.e-6) bz = 1.e-6;
	blpx = sqrt(bx * bx + by * by + bz * bz); 
	nsux = (-hfsx - x2) / bx;
	psux = ( hfsx - x2) / bx;
	nsuy = (-hfsy - y2) / by;
	psuy = ( hfsy - y2) / by;
	nsuz = (-hfsz - z2) / bz;
	psuz = ( hfsz - z2) / bz;
	dminx = (nsux < psux) ? nsux : psux;
	dminy = (nsuy < psuy) ? nsuy : psuy;
	dminz = (nsuz < psuz) ? nsuz : psuz;
	t1 = (0.0 > dminx) ? 0.0 : dminx;
	t2 = (dminy > dminz) ? dminy : dminz;
	dmin = (t1 > t2) ? t1 : t2;
	dmaxx = (nsux > psux) ? nsux : psux;
	dmaxy = (nsuy > psuy) ? nsuy : psuy;
	dmaxz = (nsuz > psuz) ? nsuz : psuz;
	t1 = (1.0 < dmaxx) ? 1.0 : dmaxx;
	t2 = (dmaxy < dmaxz) ? dmaxy : dmaxz;
	dmax = (t1 < t2) ? t1 : t2;


	if (dmax <= dmin)
	{ //fprintf(stderr, "siddon: the beam does not intersect the array (!)\n"); 
		return; }

	if (bx >= 0.0)
	{
		bilstartx = (fabs(dmin - dminx) < 1.e-6) ? 
			(int)(hfsx + dmin * bx + x2 + 0.5) : 
			(int)(hfsx + dmin * bx + x2 + 1);
		bilstopx  = (fabs(dmax - dmaxx) < 1.e-6) ? 
			(int)(hfsx + x2 + dmax * bx + 0.5) : 
			(int)(hfsx + x2 + dmax * bx);
	}
	else
	{
		bilstartx = (fabs(dmax - dmaxx) < 1.e-6) ? 
			(int)(hfsx + dmax * bx + x2 + 0.5) : 
			(int)(hfsx + dmax * bx + x2 + 1);
		bilstopx  = (fabs(dmin - dminx) < 1.e-6) ? 
			(int)(hfsx + x2 + dmin * bx + 0.5) : 
			(int)(hfsx + x2 + dmin * bx);
	}
	if (by >= 0.0)
	{
		bilstarty = (fabs(dmin - dminy) < 1.e-6) ? 
			(int)(hfsy + dmin * by + y2 + 0.5) : 
			(int)(hfsy + dmin * by + y2 + 1);
		bilstopy  = (fabs(dmax - dmaxy) < 1.e-6) ? 
			(int)(hfsy + y2 + dmax * by + 0.5) : 
			(int)(hfsy + y2 + dmax * by);
	}
	else
	{
		bilstarty = (fabs(dmax - dmaxy) < 1.e-6) ? 
			(int)(hfsy + dmax * by + y2 + 0.5) : 
			(int)(hfsy + dmax * by + y2 + 1);
		bilstopy  = (fabs(dmin - dminy) < 1.e-6) ? 
			(int)(hfsy + y2 + dmin * by + 0.5) : 
			(int)(hfsy + y2 + dmin * by);
	}
	if (bz >= 0.0)
	{
		bilstartz = (fabs(dmin - dminz) < 1.e-6) ? 
                        (int)(hfsz + dmin * bz + z2 + 0.5) : 
                        (int)(hfsz + dmin * bz + z2 + 1);
		bilstopz  = (fabs(dmax - dmaxz) < 1.e-6) ? 
                        (int)(hfsz + z2 + dmax * bz + 0.5) : 
                        (int)(hfsz + z2 + dmax * bz);
	}
	else
	{
		bilstartz = (fabs(dmax - dmaxz) < 1.e-6) ? 
                        (int)(hfsz + dmax * bz + z2 + 0.5) : 
                        (int)(hfsz + dmax * bz + z2 + 1);
		bilstopz  = (fabs(dmin - dminz) < 1.e-6) ? 
                        (int)(hfsz + z2 + dmin * bz + 0.5) : 
                        (int)(hfsz + z2 + dmin * bz);
	}

	ax[bilstartx] = (bilstartx - hfsx - x2) /bx;
	for (x = bilstartx + 1, fact = 1.0 / bx; x <= bilstopx; x++)
		ax[x] = ax[x-1] + fact;
	ay[bilstarty] = (bilstarty - hfsy - y2) /by;
	for (y = bilstarty + 1, fact = 1.0 / by; y <= bilstopy; y++)
		ay[y] = ay[y-1] + fact;
	az[bilstartz] = (bilstartz - hfsz - z2)/bz;
	for (z = bilstartz + 1, fact = 1.0 / bz; z <= bilstopz; z++)
		az[z] = az[z-1] + fact;


	if (bx > 0.0) xflag =  1, x = bilstartx;
	else          xflag = -1, x = bilstopx;
	tx = (bilstartx <= bilstopx) ? ax[x] : (ax[(x = bilstartx)] = 10.0);
	if (by > 0.0) yflag =  1, y = bilstarty;
	else          yflag = -1, y = bilstopy;
	ty = (bilstarty <= bilstopy) ? ay[y] : (ay[(y = bilstarty)] = 10.0);
	if (bz > 0.0) zflag =  1, z = bilstartz;
	else          zflag = -1, z = bilstopz;
	tz = (bilstartz <= bilstopz) ? az[z] : (az[(z = bilstartz)] = 10.0);



	BV = (bilstopx - bilstartx + 1) + (bilstopy - bilstarty + 1) + 
		(bilstopz - bilstartz + 1);
	alpha = (float *)malloc((BV+2) * sizeof(float));
	for (bv = 0; bv < BV; bv++)
	{
		if (tx <= ty && tx <= tz)
		{
			alpha[bv] = tx;
			x += xflag;
			tx = (x < bilstartx || x > bilstopx) ? 1. + dmax : ax[x];
		}
		else if (ty <= tx && ty <= tz)
		{
			alpha[bv] = ty;
			y += yflag;
			ty = (y < bilstarty || y > bilstopy) ? 1. + dmax : ay[y];
		}
		else
		{
			alpha[bv] = tz;
			z += zflag;
			tz = (z < bilstartz || z > bilstopz) ? 1. + dmax : az[z];
		}
	}

	bvattew = (float *)calloc(BV+1, sizeof(float));
	x_ = (int *)malloc((BV+1) * sizeof(int));
	y_ = (int *)malloc((BV+1) * sizeof(int));
	z_ = (int *)malloc((BV+1) * sizeof(int));

	for(bv = 1, batte = 0.0; bv < BV; bv++)
	{
		amid = (alpha[bv] + alpha[bv-1]) / 2.0;
		x_[bv] = x = (int)(hfsx + x2 + amid * bx);
		y_[bv] = y = (int)(hfsy + y2 + amid * by);
		z_[bv] = z = (int)(hfsz + z2 + amid * bz);

		if (z<0) z=0;
		if (y<0) y=0;
		if (x<0) x=0;
		if (z>=ZDim) z=ZDim-1;
		if (y>=YDim) y=YDim-1;
		if (x>=XDim) x=XDim-1;
		z_[bv]=z;
		y_[bv]=y;
		x_[bv]=x;

		if (image[x+y*XDim+z*XDim*YDim] > 1.e-15) 
		{
			blivpx = (alpha[bv] - alpha[bv-1]) * blpx; 
			bvattew_ = bvattew[bv] = blivpx;
			if (bvattew_ < 0.0)
				fprintf(stderr, "siddon: probably dimensions of ax, ay, az vectors to short !\n");
    
			gammatogammaforw += bvattew_ * image[x+y*XDim+z*XDim*YDim];
			randomforw += bvattew_ * randImage[x+y*XDim+z*XDim*YDim];
		}
	}
	if (gammatogammaforw > 0.0)
	{
		randnorm = (gammatogammaforw) / (gammatogammaforw + randomforw);
		if (randnorm > rand_cut)
		{
			gammatogammaforw = randnorm * norm / gammatogammaforw;
			for (bv = 1; bv < BV; bv++)
			{
				diffImage[x_[bv]+y_[bv]*XDim+z_[bv]*XDim*YDim] += gammatogammaforw * bvattew[bv];
				//totalImage[x_[bv]+y_[bv]*XDim+z_[bv]*XDim*YDim] += bvattew[bv];
			}
		}
	}

	free(alpha);
	free(bvattew);
	free(x_);
	free(y_);
	free(z_);
}   

void zero_out(struct image_params *ip, float *image)
{
	float PIXELLENGTH = ip->pixelLength;
	int XDim = ip->xDim;
	int YDim = ip->yDim;
	int ZDim = ip->zDim;
	
	std::set<float>	 uniqueAngles;
	for(std::multiset<float>::iterator i = ip->angles.begin(); i != ip->angles.end(); i++)
		uniqueAngles.insert(*i);

	for(std::set<float>::iterator i = uniqueAngles.begin(); i != uniqueAngles.end(); i++) {
		float theta = *i;
		float sinTheta = sinf(-theta);
		float cosTheta = cosf(-theta);

		for(int xi = 0; xi < XDim; xi++) {
		 	for(int yi = 0; yi < YDim; yi++) {
				for(int zi = 0; zi < ZDim; zi++) {
					float x = (xi - XDim/2) * PIXELLENGTH;
					float y = (yi - YDim/2) * PIXELLENGTH;
					float z = (zi - ZDim/2) * PIXELLENGTH;

					float x_rot = x;
					float y_rot = +y * cosTheta + z * sinTheta;
					float z_rot = -y * sinTheta + z * cosTheta;
					
					if ((fabs(y_rot) > ip->transaxialSize/2) || (fabs(z_rot) > ip->transaxialSize/2))
						image[xi+yi*XDim+zi*XDim*YDim] = 0;

					if(fabs(x) > ip->axialSize/2) {
						image[xi+yi*XDim+zi*XDim*YDim] = 0;
					}
				}
			}
		}
	}

/*	float radius = (ip->platesDistance - 60)/2;	
	
	for(int xi = 0; xi < XDim; xi++)
	 	for(int yi = 0; yi < YDim; yi++)
			for(int zi = 0; zi < ZDim; zi++) {
				float x = (xi - XDim/2) * PIXELLENGTH;
				float y = (yi - YDim/2) * PIXELLENGTH;
				float z = (zi - ZDim/2) * PIXELLENGTH;

				if ((y*y + z*z) > (radius*radius)) {
					image[xi+yi*XDim+zi*XDim*YDim] = 0;
				}

			}*/
}

static float gaussian(float sigma, int n, float d2) 
{
	return 1/powf(sqrtf(2*M_PI)*sigma, n) * expf(- d2/(2*sigma*sigma));
}


void gaussian_filter_3d(float *recImage, float sigma, struct image_params *ip)
{
	float PIXELLENGTH = ip->pixelLength;
	int XDim = ip->xDim;
	int YDim = ip->yDim;
	int ZDim = ip->zDim;

	float *	tempFilterResult = new float[XDim * YDim * ZDim];

	const int ks = 5;
	const int kw = 2*ks + 1;
	float weights[kw][kw][kw];
	
	for(int i = 0; i < kw; i++) {		
		for(int j = 0; j < kw; j++) {
			for (int k = 0; k < kw; k++) {
				float x = (i - ks) * PIXELLENGTH;
				float y = (j - ks) * PIXELLENGTH;
				float z = (k - ks) * PIXELLENGTH;
				float w = gaussian(sigma, 3, x*x + y*y + z*z);
				weights[i][j][k] = w;

			}
		}
	}	

	float srcTotal = 0;
	float fltTotal = 0;
	for (int xi = 0; xi < XDim; xi++) {
		for (int yi = 0; yi < YDim; yi++) {
			for (int zi = 0; zi < ZDim; zi++) {
				srcTotal += recImage[xi+yi*XDim+zi*XDim*YDim];
				float sum = 0;
				for (int k = 0; k < kw; k++) {
					if ((zi+k-ks < 0) || (zi+k-ks >= ZDim)) continue;
					for(int j = 0; j < kw; j++) {
						if ((yi+j-ks < 0) || (yi+j-ks >= YDim)) continue;
						for(int i = 0; i < kw; i++) {
							if ((xi+i-ks < 0) || (xi+i-ks >= XDim)) continue;
							sum += weights[i][j][k] * recImage[(xi+i-ks)+(yi+j-ks)*XDim+(zi+k-ks)*XDim*YDim];
						}
					}
				}
				
				tempFilterResult[xi+yi*XDim+zi*XDim*YDim] = sum;	
				fltTotal += sum;
			}
		}
	}
	
	/* 
	 * Copy and normalize 
	 */ 
	for(int i = 0; i < XDim * YDim * ZDim; i++)
		recImage[i] = tempFilterResult[i] * (srcTotal / fltTotal);
	delete [] tempFilterResult;
	
}

void gaussian_filter_21d(float *recImage, float sigma, struct image_params *ip)
{
	float PIXELLENGTH = ip->pixelLength;
	int XDim = ip->xDim;
	int YDim = ip->yDim;
	int ZDim = ip->zDim;

	float *	tempFilterResult = new float[XDim * YDim * ZDim];

	const int ks = 5;
	const int kw = 2*ks + 1;
	float weightsYZ[kw][kw];
	
	for(int j = 0; j < kw; j++) {
		for (int k = 0; k < kw; k++) {
			float y = (j - ks) * PIXELLENGTH;
			float z = (k - ks) * PIXELLENGTH;
			float w = gaussian(sigma, 2, y*y + z*z);
			weightsYZ[j][k] = w;
		}
	}	


	/*
	 * Filter Y and Z
	 */
	float srcTotal = 0;
	
	for (int xi = 0; xi < XDim; xi++) {
		for (int yi = 0; yi < YDim; yi++) {
			for (int zi = 0; zi < ZDim; zi++) {
				srcTotal += recImage[xi+yi*XDim+zi*XDim*YDim];
				float sum = 0;
				for (int k = 0; k < kw; k++) {
					if ((zi+k-ks < 0) || (zi+k-ks >= ZDim)) continue;
					for(int j = 0; j < kw; j++) {
						if ((yi+j-ks < 0) || (yi+j-ks >= YDim)) continue;
						sum += weightsYZ[j][k] * recImage[(xi)+(yi+j-ks)*XDim+(zi+k-ks)*XDim*YDim];
					}
				}				
				tempFilterResult[xi+yi*XDim+zi*XDim*YDim] = sum;	
			}
		}
	}

	/*
	 * Filter X
	 */

	float weightsX[kw];
	for(int i = 0; i < kw; i++) {
		float x = (i - ks) * PIXELLENGTH;
		float w = gaussian(sigma, 1, x*x);
		weightsX[i] = w;
	}

	float fltTotal = 0;
	for (int xi = 0; xi < XDim; xi++) {
		for (int yi = 0; yi < YDim; yi++) {
			for (int zi = 0; zi < ZDim; zi++) {				
				float sum = 0;
				for (int i = 0; i < kw; i++) {
					if ((xi+i-ks < 0) || (xi+i-ks >= XDim)) continue;
					sum += weightsX[i] * tempFilterResult[(xi+i-ks)+(yi)*XDim+(zi)*XDim*YDim];
				}							
				recImage[xi+yi*XDim+zi*XDim*YDim] = sum;	
				fltTotal += sum;
			
			}
		}
	}

	/* 
	 * Normalize 
	 */ 
	for(int i = 0; i < XDim * YDim * ZDim; i++)
		recImage[i] = recImage[i] * (srcTotal / fltTotal);
	delete [] tempFilterResult;
	
}

void gaussian_filter_2d(float *recImage, float sigma, struct image_params *ip)
{
	float PIXELLENGTH = ip->pixelLength;
	int XDim = ip->xDim;
	int YDim = ip->yDim;
	int ZDim = ip->zDim;

	float *	tempFilterResult = new float[XDim * YDim * ZDim];

	const int ks = 5;
	const int kw = 2*ks + 1;
	float weightsYZ[kw][kw];
	
	for(int j = 0; j < kw; j++) {
		for (int k = 0; k < kw; k++) {
			float y = (j - ks) * PIXELLENGTH;
			float z = (k - ks) * PIXELLENGTH;
			float w = gaussian(sigma, 2, y*y + z*z);
			weightsYZ[j][k] = w;
		}
	}	


	/*
	 * Filter Y and Z
	 */
	float srcTotal = 0;
	float fltTotal = 0;
	
	for (int xi = 0; xi < XDim; xi++) {
		for (int yi = 0; yi < YDim; yi++) {
			for (int zi = 0; zi < ZDim; zi++) {
				srcTotal += recImage[xi+yi*XDim+zi*XDim*YDim];
				float sum = 0;
				for (int k = 0; k < kw; k++) {
					if ((zi+k-ks < 0) || (zi+k-ks >= ZDim)) continue;
					for(int j = 0; j < kw; j++) {
						if ((yi+j-ks < 0) || (yi+j-ks >= YDim)) continue;
						sum += weightsYZ[j][k] * recImage[(xi)+(yi+j-ks)*XDim+(zi+k-ks)*XDim*YDim];
					}
				}				
				tempFilterResult[xi+yi*XDim+zi*XDim*YDim] = sum;	
				fltTotal += sum;
			}
		}
	}

	/* 
	 * Normalize 
	 */ 
	for(int i = 0; i < XDim * YDim * ZDim; i++)
		recImage[i] = tempFilterResult[i] * (srcTotal / fltTotal);
	delete [] tempFilterResult;
	
}
