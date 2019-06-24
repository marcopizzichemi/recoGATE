#ifndef __LMREC_METZ__DEFINED__
#define __LMREC_METZ__DEFINED__

void Metz2D(float * Data, int Xdim, int Ydim, float StdDev, float Power);
void Metz3D(float * Data, int Xdim, int Ydim, int Zdim, float StdDev, float Power);

static void metz_filter_3d(float *recImage, float sigma, float power, struct image_params *ip)
{
	
	Metz3D(recImage, ip->xDim, ip->yDim, ip->zDim, sigma, power);

}


#endif
