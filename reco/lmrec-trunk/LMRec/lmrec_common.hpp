#ifndef __LMREC_COMMON_HPP__DEFINED__
#define __LMREC_COMMON_HPP__DEFINED__

#include <vector>
#include <set>
#include "../lmFormats/formats.hpp"

#include <boost/random.hpp>
#include <boost/nondet_random.hpp>

struct image_params {
	float transaxialSize;
	float axialSize;
	float platesDistance;
	float pixelLength;
	std::multiset<float> angles;
	int xDim;
	int yDim;
	int zDim;
};

struct rec_kernel_parameter {
	struct image_params *ip;

	DKFZFormat* listData;	
	int beginNum;
	int totalNum;	
	float* diffImage;
	float* normImage1;
	float* normImage2;
	float* recImage;
	float* randImage;
	float rand_cut;

	int nRays;
	boost::mt19937 *generator;
};	


void* rec_kernel(void* parameter);

void siddon_rec1(struct image_params *ip, 
	float x1, float y1, float z1, float x2, float y2, float z2, 
	float *image, float *diffImage, float norm, float *randImage, float rand_cut);

void zero_out(struct image_params *ip, float *image);

void gaussian_filter_3d(float *recImage, float sigma, struct image_params *ip);
void gaussian_filter_21d(float *recImage, float sigma, struct image_params *ip);
void gaussian_filter_2d(float *recImage, float sigma, struct image_params *ip);


#endif 
