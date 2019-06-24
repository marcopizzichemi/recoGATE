#ifndef __PEMGEOMETRY_HPP__DEFINED__
#define __PEMGEOMETRY_HPP__DEFINED__

#include <math.h>

static const float PEM_angle_tolerance = 0.5 * M_PI/180;

static const unsigned PEM_uDim = 172;
static const float PEM_xLength = 171.5 + 2.0;
static const unsigned PEM_nCrystals_x = 64;
static const float PEM_crystals_x[PEM_nCrystals_x] = {
	
	0, 2.3, 4.6, 6.9, 9.2, 11.5, 13.8, 16.1,
	21.9, 24.2, 26.5, 28.8, 31.1, 33.4, 35.7, 38,
	44.5, 46.8, 49.1, 51.4, 53.7, 56, 58.3, 60.6,
	66.4, 68.7, 71, 73.3, 75.6, 77.9, 80.2, 82.5,
	89, 91.3, 93.6, 95.9, 98.2, 100.5, 102.8, 105.1,
	110.9, 113.2, 115.5, 117.8, 120.1, 122.4, 124.7, 127,
	133.5, 135.8, 138.1, 140.4, 142.7, 145, 147.3, 149.6,
	155.4, 157.7, 160, 162.3, 164.6, 166.9, 169.2, 171.5
};

static const unsigned PEM_vDim = 151;
static const float PEM_yLength = 150.2 + 2.0;
static const unsigned PEM_nCrystals_y = 48;
static const float PEM_crystals_y[PEM_nCrystals_y] = {
	0, 2.3, 4.9, 7.2, 13, 15.3, 17.9, 20.2,
	26, 28.3, 30.9, 33.2, 39, 41.3, 43.9, 46.2,
	52, 54.3, 56.9, 59.2, 65, 67.3, 69.9, 72.2,
	78, 80.3, 82.9, 85.2, 91, 93.3, 95.9, 98.2,
	104, 106.3, 108.9, 111.2, 117, 119.3, 121.9, 124.2,
	130, 132.3, 134.9, 137.2, 143, 145.3, 147.9, 150.2
};
	
#endif
