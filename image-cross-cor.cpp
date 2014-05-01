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

g++ -lgdal1.6.0 -lproj -lfftw3 -lfftw3_threads -lm image.cpp -o image -I/usr/include/gdal
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

void error(const char* msg)
{
	throw(std::out_of_range(msg));
	exit(0);
} 

#define MY_X(ii, jj) x[((ii)*cols) + (jj)]
void fft2shift(unsigned char *x, int rows, int cols )
{
	int m2, n2;
	int i, k;
	unsigned char tmp1, tmp2;

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

inline void runFFT(fftwf_plan plan, GDALDataset *srcDS, complex float *img, int band, GDALDataset *dstDS)
{
	const size_t px_count = srcDS->GetRasterXSize() * srcDS->GetRasterYSize();

	/* Note: sizeof(complex float) * dstDS->GetRasterXSize() is the length of a scanline
         *       in the destination image's buffer. This is used in case the source is smaller
         *       than the destination (as we don't want to rescale... I think.
    	 */
	srcDS->GetRasterBand(band)->RasterIO( GF_Read, 0, 0, 
				   srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
				   img, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
				   GDT_CFloat32, 0, sizeof(complex float) * dstDS->GetRasterXSize());
	fftwf_execute(plan);

	complex float norm = csqrt(px_count + 0I);
	for(int i = 0; i < px_count; i++) {
		img[i] = img[i] / norm;
	}
}


float* CalcFFT(GDALDataset *srcDS1, GDALDataset *srcDS2, GDALDataset *dstDS) 
{
	fftwf_plan plan1, plan2, planI;
	fftwf_complex *img1, *img2;
	unsigned char *out;
	int band;

	const size_t px_count = dstDS->GetRasterXSize() * dstDS->GetRasterYSize();
	const size_t buffer_len = sizeof(fftwf_complex) * px_count;
	img1  = (fftwf_complex*) fftwf_malloc(buffer_len);
	img2  = (fftwf_complex*) fftwf_malloc(buffer_len);
	out   = (unsigned char*) fftwf_malloc(sizeof(unsigned char) * px_count); 
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
		runFFT( plan1, srcDS1, img1, band, dstDS );
		printf("FFT 2 band %d\n", band);
		runFFT( plan2, srcDS2, img2, band, dstDS );

		
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
		/* img2 should now be real - normalize 0-255 and -- write output */
		printf("Save band %d; min = %f max = %f\n", band, min, max);
		for(int i = 0; i < px_count; i++) {
			out[i] = floor( ((cabs(img2[i]) - min) / (max-min) ) * 255.0 );
		}

		fft2shift(out, dstDS->GetRasterYSize(), dstDS->GetRasterXSize());

 		dstDS->GetRasterBand(band)->RasterIO( GF_Write, 0, 0, 
					   dstDS->GetRasterXSize(), dstDS->GetRasterYSize(),
					   out, dstDS->GetRasterXSize(), dstDS->GetRasterYSize(),
					   GDT_Byte, 0, 0);
	}



	fftwf_destroy_plan(plan1);
	fftwf_destroy_plan(plan2);
	fftwf_destroy_plan(planI);
	fftwf_free(img1);
	fftwf_free(img2);
	fftwf_free(out);
}

int ConvertToPolar()
{

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
	
	GDALDataset *dstDS = gd->Create(fileName, width, height, bands, GDT_Byte, options);
	
	return dstDS;
}


int main(int argc, char** argv)
{
	char* srcFileName1;
	char* srcFileName2;
	char* dstFileName;

	GDALDataset* srcDS1, *srcDS2;
	GDALDataset* dstDS;

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

	int x_size, y_size, n_bands;
	x_size = MY_MAX(srcDS1->GetRasterXSize(), srcDS2->GetRasterXSize());
	y_size = MY_MAX(srcDS1->GetRasterYSize(), srcDS2->GetRasterYSize());
	n_bands = MY_MIN(srcDS1->GetRasterCount(), srcDS2->GetRasterCount());

	dstDS = CreateOutputDataset(dstFileName, 
	                            x_size, y_size, n_bands
	                           );

	if(dstDS == NULL)
		error("Could not create destination dataset");

	CalcFFT(srcDS1, srcDS2, dstDS);

	delete srcDS1;
	delete srcDS2;
	delete dstDS;

	return 0;	
}
