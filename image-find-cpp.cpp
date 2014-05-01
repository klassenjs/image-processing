/*
 * Finds the Conjucate Principle Point on an image
 * made with image-cross-cor

 * Copyright (c) 2011-2013, James Klassen
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <stdexcept>
#include <string>
#include <proj_api.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "/usr/include/complex.h"
#include <math.h>
#include <fftw3.h>

void error(const char* msg)
{
	throw(std::out_of_range(msg));
	exit(0);
} 

void CreateOutputDataset(char* fileName, int width, int height, float *img)
{
	GDALDriverManager *gdm = GetGDALDriverManager();
	if(gdm == NULL)
		error("GetGDALDriverManager() failed!");

	GDALDriver *gd = gdm->GetDriverByName("GTiff");
	if(gd == NULL)
		error("Get GTiff Driver failed!");

	char* options[2];
	options[0] = "INTERLEAVE=BAND";
	options[1] = NULL;
	
	GDALDataset *dstDS = gd->Create(fileName, width, height, 1, GDT_Float32, options);
 	dstDS->GetRasterBand(1)->RasterIO( GF_Write, 0, 0, 
					   dstDS->GetRasterXSize(), dstDS->GetRasterYSize(),
					   img, dstDS->GetRasterXSize(), dstDS->GetRasterYSize(),
					   GDT_Float32, 0, 0);


	delete dstDS;	
}

void findCPPMax(GDALDataset *srcDS, float mask_pct) 
{
	float *img, *buf;

	const int rows = srcDS->GetRasterYSize();
	const int cols = srcDS->GetRasterXSize();
	const int bands = srcDS->GetRasterCount();

	const size_t px_count = rows * cols;
	img = (float*) fftwf_malloc(sizeof(float) * px_count); 
	buf = (float*) fftwf_malloc(sizeof(float) * px_count); 
		/* ^ not used in fft, but aligned is good anyway */
	if(img == NULL || buf == NULL)
		error("Could not allocate memory\n");

	if(srcDS->GetRasterCount() < 1)
		error("Image needs to have at least one band\n");


	/* Multiply all bands together into composite image */
 	srcDS->GetRasterBand(1)->RasterIO( GF_Read, 0, 0, 
			   cols, rows, img, cols, rows, GDT_Float32, 0, 0);

	for(int band = 2; band <= bands; band++) {
 		srcDS->GetRasterBand(band)->RasterIO( GF_Read, 0, 0, 
				   cols, rows, buf, cols, rows, GDT_Float32, 0, 0);

		for(int i = 0; i < px_count; i++) {
			img[i] = img[i] * buf[i];
		}
	}

	/* Row and column of image center */
	const int pp_row = rows / 2;
	const int pp_col = cols / 2;

	const int rows2 = pp_row * 2; /* Avoid odd pixels at end */
	const int cols2 = pp_col * 2; /* image-cross-cor fft2shift has a bug with this */

	/* Calculate Mask Rectangle */
	if(mask_pct > 99) mask_pct = 99;
	if(mask_pct <= 1) mask_pct = 0;

	const int mask_min_row = pp_row - ( pp_row * mask_pct / 100);
	const int mask_max_row = pp_row + ( pp_row * mask_pct / 100);

	const int mask_min_col = pp_col - ( pp_col * mask_pct / 100);
	const int mask_max_col = pp_col + ( pp_col * mask_pct / 100);

	/* Find peak response that is outside the mask */
	float max = -1;
	int max_col = 0, max_row = 0;

	for(int row = 0; row < rows2; row++) {
		int row_offset = row * cols;

		if((row >= mask_min_row) && (row < mask_max_row)) {  /* we need to skip central columns */
			for(int col = 0; col < mask_min_col; col++) {
				float v = img[row_offset + col];
				if(v > max) {
					max = v;
					max_col = col;
					max_row = row;
				}
			}
			for(int col = mask_max_col; col < cols2; col++) {
				float v = img[row_offset + col];
				if(v > max) {
					max = v;
					max_col = col;
					max_row = row;
				}
			}
		} else {
			for(int col = 0; col < cols2; col++) {
				float v = img[row_offset + col];
				if(v > max) {
					max = v;
					max_col = col;
					max_row = row;
				}
			}
		}
	}

	//CreateOutputDataset("out-cpp.tif", cols, rows, img); /* debug so we can see what it's doing */
	printf("\t%d\t%d\t%f\n", max_col - pp_col, max_row - pp_row, max);

	fftwf_free(img);
	fftwf_free(buf);
}

int main(int argc, char** argv)
{
	char* srcFileName;
	float mask_pct;

	GDALDataset* srcDS;

	if(argc != 3)
		error("Usage: prog infile mask_percent ");

	srcFileName = argv[1];
	mask_pct = atof(argv[2]);

	GDALAllRegister();

	srcDS = (GDALDataset*) GDALOpen( srcFileName, GA_ReadOnly );

	if(srcDS == NULL)
		error("Could not open source dataset");

	//printf("Got image %d x %d x %d\n",
	//      srcDS->GetRasterXSize(),
	//      srcDS->GetRasterYSize(),
	//      srcDS->GetRasterCount()
	//);


	findCPPMax(srcDS, mask_pct);

	delete srcDS;

	return 0;	
}
