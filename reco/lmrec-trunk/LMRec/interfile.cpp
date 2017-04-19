#include "interfile.hpp"
#include <assert.h>
#include <stdio.h>

void writeInterFile(char *prefix, struct image_params *ip, float *image)
{

	char fName[1024];
	sprintf(fName, "%s.hv", prefix);
	FILE * header = fopen(fName, "w");
	assert(header != NULL);

	sprintf(fName, "%s.v", prefix);
	FILE * data = fopen(fName, "wb");
	assert(data != NULL);

	int XDim = ip->xDim;
	int YDim = ip->yDim;
	int ZDim = ip->zDim;


	fprintf(header, "! INTERFILE  :=\n");
	fprintf(header, "name of data file := %s\n", fName);
	fprintf(header, "!GENERAL DATA :=\n");
	fprintf(header, "!GENERAL IMAGE DATA :=\n");
	fprintf(header, "!type of data := PET\n");
	fprintf(header, "imagedata byte order := LITTLEENDIAN\n");
	fprintf(header, "!PET STUDY (General) :=\n");
	fprintf(header, "!PET data type := Image\n");
	fprintf(header, "process status := Reconstructed\n");
	fprintf(header, "!number format := float\n");
	fprintf(header, "!number of bytes per pixel := 4\n");
	fprintf(header, "number of dimensions := 3\n");
	fprintf(header, "matrix axis label [1] := x\n");
	fprintf(header, "!matrix size [1] := %d\n", ZDim);
	fprintf(header, "scaling factor (mm/pixel) [1] := %f\n", ip->pixelLength);
	fprintf(header, "matrix axis label [2] := y\n");
	fprintf(header, "!matrix size [2] := %d\n", YDim);
	fprintf(header, "scaling factor (mm/pixel) [2] := %f\n", ip->pixelLength);
	fprintf(header, "matrix axis label [3] := z\n");
	fprintf(header, "!matrix size [3] := %d\n", XDim);
	fprintf(header, "scaling factor (mm/pixel) [3] := %f\n", ip->pixelLength);
	fprintf(header, "number of time frames := 1\n");
	fprintf(header, "image scaling factor[1] := 1\n");
	fprintf(header, "data offset in bytes[1] := 0\n");
	fprintf(header, "quantification units := 1\n");
	fprintf(header, "!END OF INTERFILE :=\n");

	for(int xi = 0; xi < XDim; xi++)
		for(int zi = 0; zi < ZDim; zi++)
			for(int yi = 0; yi < YDim; yi++)
				fwrite((void*)&image[(XDim-xi-1)+(ZDim-zi-1)*XDim+(yi)*XDim*YDim], sizeof(float), 1, data);



	fclose(data);
	fclose(header);


}
