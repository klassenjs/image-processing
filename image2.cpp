/*
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

g++ -I/Applications-POSIX/osgeo/include -I/opt/local/include -L/Applications-POSIX/osgeo/lib -L/opt/local/lib -lgdal -lproj -lfftw3f -lfftw3f_threads -lm -O2 image.cpp -o image

g++ -lgdal1.6.0 -lproj -lfftw3 -lfftw3_threads -lm image2.cpp -o image2 -I/usr/include/gdal
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

#ifndef CORES
 #define CORES 2
#endif
#define MY_MIN(a,b) (a<b)? a : b
#define MY_MAX(a,b) (a>b)? a : b


struct Rect {
	int minx;
	int miny;
	int width;
	int height;
};


void error(const char* msg)
{
	throw(std::out_of_range(msg));
	exit(0);
} 

#define MY_X(ii, jj) x[((ii)*cols) + (jj)]
void fft2shift(float *x, int rows, int cols )
{
	int m2, n2;
	int i, k;
	float tmp1, tmp2;

	m2 = rows / 2;    // half of row dimension
	n2 = cols / 2;    // half of column dimension

	for (i = 0; i < rows; i++)  // row
	{
		for(k = 0; k < n2; k++) // col
		{
			tmp1 = MY_X(i, k + n2);
			MY_X(i, k + n2) = MY_X(i, k);
			MY_X(i, k) = tmp1;
		}
	}
	for (i = 0; i < m2; i++)  /* row */
	{
		for(k = 0; k < cols; k++) /* col */
		{
			tmp2 = MY_X(i + m2, k);
			MY_X(i + m2, k) = MY_X(i, k);
			MY_X(i, k) = tmp2;
		}
	}

}
#undef MY_X

inline void runFFT(fftwf_plan plan, GDALDataset *srcDS, complex float *img, int band, GDALDataset *dstDS, struct Rect bbox)
{
	const size_t px_count = dstDS->GetRasterXSize() * dstDS->GetRasterYSize();

	/* Note: sizeof(complex float) * dstDS->GetRasterXSize() is the length of a scanline
         *       in the destination image's buffer. This is used in case the source is smaller
         *       than the destination (as we don't want to rescale... I think.)
    	 */
	srcDS->GetRasterBand(band)->RasterIO( GF_Read, bbox.minx, bbox.miny, 
				   bbox.width, bbox.height,
				   img, bbox.width, bbox.height,
				   GDT_CFloat32, 0, sizeof(complex float) * dstDS->GetRasterXSize());
	fftwf_execute(plan);

	complex float norm = csqrt(px_count + 0I);
	for(int i = 0; i < px_count; i++) {
		img[i] = img[i] / norm;
	}
}


float* CalcFFT(GDALDataset *srcDS1, GDALDataset *srcDS2, GDALDataset *dstDS, struct Rect bbox1, struct Rect bbox2) 
{
	fftwf_plan plan1, plan2, planI;
	fftwf_complex *img1, *img2;
	float *out;
	int band;

	const size_t px_count = dstDS->GetRasterXSize() * dstDS->GetRasterYSize();
	const size_t buffer_len = sizeof(fftwf_complex) * px_count;
	img1  = (fftwf_complex*) fftwf_malloc(buffer_len);
	img2  = (fftwf_complex*) fftwf_malloc(buffer_len);
	out   = (float*) fftwf_malloc(sizeof(float) * px_count); 
		/* ^ not used in fft, but aligned is good anyway */
	if(img1 == NULL || img2 == NULL || out == NULL)
		error("Could not allocate memory\n");

	if(fftwf_init_threads())
		fftwf_plan_with_nthreads(CORES);

	plan1 = fftwf_plan_dft_2d(dstDS->GetRasterYSize(), dstDS->GetRasterXSize(), 
	                        img1, img1, FFTW_FORWARD, FFTW_ESTIMATE);

	plan2 = fftwf_plan_dft_2d(dstDS->GetRasterYSize(), dstDS->GetRasterXSize(), 
	                        img2, img2, FFTW_FORWARD, FFTW_ESTIMATE);

	planI = fftwf_plan_dft_2d(dstDS->GetRasterYSize(), dstDS->GetRasterXSize(), 
	                        img2, img2, FFTW_BACKWARD, FFTW_ESTIMATE);

	if(plan1 == NULL || plan2 == NULL || planI == NULL)
		error("Could not plan FFT\n");

	for(band = 1; band <= dstDS->GetRasterCount(); band++) {
		printf("FFT 1 band %d\n", band);
		runFFT( plan1, srcDS1, img1, band, dstDS, bbox1 );
		printf("FFT 2 band %d\n", band);
		runFFT( plan2, srcDS2, img2, band, dstDS, bbox2 );

		
		printf("Complex Conj band %d\n", band);
		/* mult img1 and conj of img2 */
		for(int px = 0; px < px_count; px++) {
			img2[px] = img1[px] * conj(img2[px]);
		}
	
		/* IFFT of result */	
		printf("IFFT band %d\n", band);
		fftwf_execute(planI);	

		printf("normalize band %d\n", band);
		complex float norm = csqrt(px_count + 0I);
		float max = cabs(img2[0] / norm);
		float min = cabs(img2[0] / norm);
		for(int i = 0; i < px_count; i++) {
			img2[i] = img2[i] / norm;
			
			if(cabs(img2[i]) < min)
				min = cabs(img2[i]);
			if(cabs(img2[i]) > max)
				max = cabs(img2[i]);
		}
		/* img2 should now be real - normalize 0.0-1.0 and -- write output */
		printf("Save band %d; min = %f max = %f\n", band, min, max);
		for(int i = 0; i < px_count; i++) {
			out[i] = ((cabs(img2[i]) - min) / (max-min) );
		}

		fft2shift(out, dstDS->GetRasterYSize(), dstDS->GetRasterXSize());

 		dstDS->GetRasterBand(band)->RasterIO( GF_Write, 0, 0, 
					   dstDS->GetRasterXSize(), dstDS->GetRasterYSize(),
					   out, dstDS->GetRasterXSize(), dstDS->GetRasterYSize(),
					   GDT_Float32, 0, 0);
	}

	fftwf_destroy_plan(plan1);
	fftwf_destroy_plan(plan2);
	fftwf_destroy_plan(planI);
	fftwf_free(img1);
	fftwf_free(img2);
	fftwf_free(out);
}


GDALDataset* CreateOutputDataset(char* fileName, int width, int height, int bands)
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
	
	GDALDataset *dstDS = gd->Create(fileName, width, height, bands, GDT_Float32, options);
	
	return dstDS;
}


float findCPPMax(GDALDataset *srcDS, float mask_pct, int* offset_x, int* offset_y) 
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

	*offset_x = max_col - pp_col;
	*offset_y = max_row - pp_row;

	fftwf_free(img);
	fftwf_free(buf);

	return max;
}


void newBBox(int quadrant, const struct Rect in, struct Rect* out)
{
	out->width = in.width / 2;
	out->height = in.height / 2;
	switch(quadrant) {
		case 1:
			out->minx = in.minx + (in.width)/2;
			out->miny = in.miny + (in.height)/2;
			break;
		case 2:
			out->minx = in.minx;
			out->miny = in.miny + (in.height)/2;
			break;
		case 3:
			out->minx = in.minx;
			out->miny = in.miny;
			break;
		case 4:
			out->minx = in.minx + (in.width)/2;
			out->miny = in.miny;
			break;
		default:
			out->minx = in.minx;
			out->miny = in.miny;
	}
}

int calcDisparityMap(GDALDataset* srcDS1, GDALDataset* srcDS2, const struct Rect bbox1, const struct Rect bbox2, const int initial_offset_x, const int initial_offset_y, const char* dstFileName)
{
	GDALDataset* dstDS;
	char* newDstFileName;

	struct Rect r1, r2;

	int offset_x = 0, offset_y = 0;
	int x_size, y_size, n_bands;

	x_size = MY_MAX(bbox1.width, bbox2.width);
	y_size = MY_MAX(bbox1.height, bbox2.height);
	n_bands = MY_MIN(srcDS1->GetRasterCount(), srcDS2->GetRasterCount());

	printf("Processing: (%d,%d) w: %d h: %d\n", bbox1.minx, bbox1.miny, bbox1.width, bbox1.height);



	newDstFileName = (char*)malloc(strlen(dstFileName) + 100);
	sprintf(newDstFileName, "%s_%d-%d_%d-%d.tif", dstFileName, bbox1.minx, bbox1.miny, bbox1.width, bbox1.height);

	dstDS = CreateOutputDataset(newDstFileName, 
	                            x_size, y_size, n_bands
	                           );
	free(newDstFileName);

	if(dstDS == NULL)
		error("Could not create destination dataset");


	CalcFFT(srcDS1, srcDS2, dstDS, bbox1, bbox2);
	float max = findCPPMax(dstDS, 0, &offset_x, &offset_y);
	delete dstDS;

	offset_x += initial_offset_x;
	offset_y += initial_offset_y;

	printf("\t%d\t%d\t%f\n", offset_x, offset_y, max);

	if(bbox1.width > 160 && bbox1.height > 160 && bbox2.width > 160 && bbox2.height > 160 && max > 0.75) {
		// Quarter and try compute again.

		for(int i = 1; i <= 4; i++) {
			newBBox(i, bbox1, &r1);	
			newBBox(i, bbox2, &r2);
			r2.minx += offset_x;
			r2.miny += offset_y;

			calcDisparityMap(srcDS1, srcDS2, r1, r2, offset_x, offset_y, dstFileName);
		}
	} else {
		// Final Answer
		printf("(%d, %d) offset (%d, %d)\n", bbox1.minx + bbox1.width/2, bbox1.miny + bbox1.width/2, offset_x, offset_y);
	}

	return 0;
}

int main(int argc, char** argv)
{
	char* srcFileName1;
	char* srcFileName2;
	char* dstFileName;

	GDALDataset* srcDS1, *srcDS2;
	
	if(argc != 4)
		error("Usage: prog infile1 infile2 outfile");

	srcFileName1 = argv[1];
	srcFileName2 = argv[2];
	dstFileName = argv[3];

	GDALAllRegister();

	srcDS1 = (GDALDataset*) GDALOpen( srcFileName1, GA_ReadOnly );
	srcDS2 = (GDALDataset*) GDALOpen( srcFileName2, GA_ReadOnly );

	if(srcDS1 == NULL || srcDS2 == NULL)
		error("Could not open source dataset");

	printf("Got image %d x %d x %d\n",
	      srcDS1->GetRasterXSize(),
	      srcDS1->GetRasterYSize(),
	      srcDS1->GetRasterCount()
	);
	printf("Got image %d x %d x %d\n",
	      srcDS2->GetRasterXSize(),
	      srcDS2->GetRasterYSize(),
	      srcDS2->GetRasterCount()
	);
	if((srcDS1->GetRasterXSize() != srcDS2->GetRasterXSize()) || 
           (srcDS1->GetRasterYSize() != srcDS2->GetRasterYSize()) ||
	   (srcDS1->GetRasterCount() != srcDS2->GetRasterCount()))
		fprintf(stderr, "Warning: Source dataset geometries should match!\n");

	struct Rect bbox1, bbox2;
	bbox1.minx = bbox1.miny = bbox2.minx = bbox2.miny = 0;
	bbox1.width = srcDS1->GetRasterXSize();
	bbox1.height = srcDS1->GetRasterYSize();
	bbox2.width = srcDS2->GetRasterXSize();
	bbox2.height = srcDS2->GetRasterYSize();

	calcDisparityMap(srcDS1, srcDS2, bbox1, bbox2, 0, 0, dstFileName);

	delete srcDS1;
	delete srcDS2;
	

	return 0;	
}
