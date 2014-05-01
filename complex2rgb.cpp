/*
 *  complex2rgb.cpp
 *  
 *
 *  Created by James Klassen on 2011/2/9.
 *
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


/*
 g++ -I/Applications-POSIX/osgeo/include -L/Applications-POSIX/osgeo/lib -lgdal -lm -O2 complex2rgb.cpp -o complex2rgb
 
 */

#include <stdexcept>
#include <string>
#include <gdal_priv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "/usr/include/complex.h"
#include <math.h>

void error(const char* msg)
{
	throw(std::out_of_range(msg));
} 


//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7
// http://alienryderflex.com/quicksort/
void quickSort(float *arr, int elements) {
	
#define  MAX_LEVELS  300
	int    i=0, L, R, beg[MAX_LEVELS], end[MAX_LEVELS], swap ;
	float  piv ;
	
	beg[0]=0; end[0]=elements;
	while (i>=0) {
		L=beg[i]; R=end[i]-1;
		if (L<R) {
			piv=arr[L];
			while (L<R) {
				while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
				while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L]; }
			arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
				swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
				swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; }}
		else {
			i--; 
		}
	}
}

inline void Complex2RGB(float mag, float phase, char* rgb)
{
	/* http://snipplr.com/view/14590/hsv-to-rgb/ */
	float r,g,b;
	
	float h = phase * (180 / 3.14159) / 60.0;
	float s = 1.0;
	float v = mag;
	
	int i = floor(h);
	float f = h - i;
	float p = v * (1 - s);
	float q = v * (1 - s * f);
	float t = v * (1 - s * (1 - f));
	switch(i) {
		case 0:
			r = v; g = t; b = p;
			break;
		case 1:
			r = q; g = v; b = p;
			break;
		case 2:
			r = p; g = v; b = t;
			break;
		case 3:
			r = p; g = q; b = v;
			break;
		case 4:
			r = t; g = p; b = v;
			break;
		default:
			r = v; g = p; b = q;
	}
	rgb[0] = floor(r * 255);
	rgb[1] = floor(g * 255);
	rgb[2] = floor(b * 255);
}


int ConvertToRGB(GDALDataset *srcDS, GDALDataset *dstDS, int band) 
{
	const size_t px_count = srcDS->GetRasterXSize() * srcDS->GetRasterYSize();
	const size_t buffer_len = sizeof(complex float) * px_count;
	
	printf("Allocating %ld MB\n", buffer_len / (1024 * 1024));
	complex float* in  = (complex float*) malloc(buffer_len);
	
	printf("Processing band %d\n", band);
	srcDS->GetRasterBand(1)->RasterIO( GF_Read, 0, 0, 
									  srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  in, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  GDT_CFloat32, 0, 0);
	

	for(int i = 0; i < 10; i++) {
		printf("%f + %fj ", creal(in[i]), cimag(in[i]));
	}
	printf("\n");
	
	printf("Normalize bands\n");
	float min = cabs(in[0]);
	float max = cabs(in[0]);
	float *mag = (float*) malloc(sizeof(float) * px_count);
	for(int px ; px < px_count; px++) {
		mag[px] = cabs(in[px]);
		if(mag[px] > max) max = mag[px];
		if(mag[px] < min) min = mag[px];
	}
	
	//printf("Call quicksort\n");
	//quickSort(mag, px_count);
	//for(int i = 0; i < px_count; i++)
	//	printf("%d\t%f\n", i, mag[i]);

	
	//printf("Calc bins\n");
	//float bins[256];
	//for(int bin = 0; bin < 256; bin++) {
	//	int i = floor((px_count-1) * ((float)bin / 255.0));
	//	bins[bin] = mag[i];
	//	
	//	printf("bin %d = %f (i = %d of %ld)\n", bin, mag[i], i, px_count);
	//}

	free(mag);

	printf("Convert colorspace\n");
	char* r_band = (char*)malloc(sizeof(char) * px_count);
	char* g_band = (char*)malloc(sizeof(char) * px_count);
	char* b_band = (char*)malloc(sizeof(char) * px_count);
	
	char rgb[3];
	for(int px = 0; px < px_count; px++) {
		float mag, phase;
		/* mag = log | DN | (dB) */
		/* mag has to be between 0 and 1 */
		mag = log10f( (10 * cabs(in[px]) ) / max );
		if(mag < 0) mag = 0;

		mag = ( cabs(in[px]) - min ) / ( max - min );
		/* Now mag is scaled 0..1 */
		mag = log10f( ((2*mag) + 1) / 0.30 );
	
		phase = carg(in[px]);
		
		//int bin = 0;
		//while(bin < 255) {
		//	if(bins[bin] > mag) {
		//		break;
		//	}
		//	bin++;
		//}
		//mag = bin / 256.0;
		
		//if(px < 100)
		//	printf("%f\n", mag);
		Complex2RGB(mag, phase, rgb);
		r_band[px] = rgb[0];
		g_band[px] = rgb[1];
		b_band[px] = rgb[2];
		//r_band[px] = (mag * 256);
		//g_band[px] = (mag * 256);
		//b_band[px] = (mag * 256);
	}
	
	dstDS->GetRasterBand(1)->RasterIO( GF_Write, 0, 0,
									  srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  r_band, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  GDT_Byte, 0, 0 );
	dstDS->GetRasterBand(2)->RasterIO( GF_Write, 0, 0,
									  srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  g_band, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  GDT_Byte, 0, 0 );
	dstDS->GetRasterBand(3)->RasterIO( GF_Write, 0, 0,
									  srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  b_band, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
									  GDT_Byte, 0, 0 );
	
	free(r_band);
	free(g_band);
	free(b_band);
	free(in);
}

int main(int argc, char** argv)
{
	char* srcFileName;
	char* dstFileName;
	
	GDALDataset* srcDS;
	GDALDataset* dstDS;
	
	if(argc != 4)
		error("Usage: prog infile outfile band");
	
	srcFileName = argv[1];
	dstFileName = argv[2];
	int band = atoi(argv[3]);
	
	GDALAllRegister();
	
	srcDS = (GDALDataset*) GDALOpen( srcFileName, GA_ReadOnly );
	
	if(srcDS == NULL)
		error("Could not open source dataset");
	
	printf("Got image %d x %d x %d\n",
		   srcDS->GetRasterXSize(),
		   srcDS->GetRasterYSize(),
		   srcDS->GetRasterCount()
		   );
	
	GDALDriverManager *gdm = GetGDALDriverManager();
	if(gdm == NULL)
		error("GetGDALDriverManager() failed!");
	
	GDALDriver *gd = gdm->GetDriverByName("GTiff");
	if(gd == NULL)
		error("Get GTiff Driver failed!");
	
	char* options[2];
	options[0] = "INTERLEAVE=BAND";
	options[1] = NULL;
	dstDS = gd->Create(dstFileName, srcDS->GetRasterXSize(),
					   srcDS->GetRasterYSize(), 3, GDT_Byte, options);
	
	if(dstDS == NULL)
		error("Could not create destination dataset");
	
	ConvertToRGB(srcDS, dstDS, band);
	
	delete srcDS;
	delete dstDS;
	
	return 0;	
}
