/*
g++ -I/Applications-POSIX/osgeo/include -I/opt/local/include -L/Applications-POSIX/osgeo/lib -L/opt/local/lib -lgdal -lproj -lfftw3f -lfftw3f_threads -lm -O2 image.cpp -o image

g++ -lgdal1.6.0 -lproj -lfftw3f -lfftw3f_threads -lm image3.cpp -o image3 -I/usr/include/gdal
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



#define MY_X(ii, jj) x[((ii)*width) + (jj)]
void fft2shift(float *x, int width, int height )
{
	int m2, n2;
	int i, k;
	float tmp1, tmp2;

	m2 = height / 2;    // half of row dimension
	n2 = width / 2;    // half of column dimension

	for (i = 0; i < height; i++)  // row
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
		for(k = 0; k < width; k++) /* col */
		{
			tmp2 = MY_X(i + m2, k);
			MY_X(i + m2, k) = MY_X(i, k);
			MY_X(i, k) = tmp2;
		}
	}

}
#undef MY_X

inline void runFFT(fftwf_plan plan, GDALDataset *srcDS, complex float *img, int band, struct Rect bbox, int width, int height)
{
	const size_t px_count = width * height;

	/* Need to pad src bbox if out of bounds */
	//int src_width = srcDS->GetRasterXSize();
	//int src_height = srcDS->GetRasterYSize();
	
	//struct Rect src_bbox;
	
	//src_bbox.minx = MY_MAX(bbox.minx, 0);
	//src_bbox.miny = MY_MAX(bbox.miny, 0);
	//int maxx = MY_MIN(bbox.minx + bbox.width, src_width);
	//int maxy = MY_MIN(bbox.miny + bbox.height, src_height);
	//src_bbox.width = maxx - src_bbox.minx;
	//src_bbox.height = maxy - src_bbox.miny;
	
	//int x_offset = 0;//MY_MAX(src_bbox.minx - bbox.minx, 0);
	//int y_offset = 0;//MY_MAX(src_bbox.miny - bbox.miny, 0);
	
	//printf("minx: %d, miny: %d, width: %d, height %d, x_offset: %d, y_offset: %d\n",
	//       src_bbox.minx, src_bbox.miny, src_bbox.width, src_bbox.height, x_offset, y_offset);

	//complex float* img2 = img + width * y_offset;

	for(int i = 0 ; i < px_count; i++)
	  img[i] = 0.0;

	srcDS->GetRasterBand(band)->RasterIO( GF_Read, bbox.minx, bbox.miny, 
				   bbox.width, bbox.height,
				   img, bbox.width, bbox.height,
				   GDT_CFloat32, 0, sizeof(complex float) * width);
	fftwf_execute(plan);

	complex float norm = csqrt(px_count + 0I);
	for(int i = 0; i < px_count; i++) {
		img[i] = img[i] / norm;
	}
}

float findOffset(GDALDataset *srcDS1, GDALDataset *srcDS2, struct Rect bbox1, struct Rect bbox2, int width, int height, int* offset_x, int* offset_y) 
{
	fftwf_plan plan1, plan2, planI;
	fftwf_complex *img1, *img2;
	float *out;
	int band;
	int n_bands = MY_MIN(srcDS1->GetRasterCount(), srcDS2->GetRasterCount());
	n_bands = 1;
	const size_t px_count = width * height;
	const size_t buffer_len = sizeof(fftwf_complex) * px_count;

	img1  = (fftwf_complex*) fftwf_malloc(buffer_len);
	img2  = (fftwf_complex*) fftwf_malloc(buffer_len);
	out   = (float*) fftwf_malloc(sizeof(float) * px_count);
		/* ^ not used in fft, but aligned is good anyway */
	if(img1 == NULL || img2 == NULL || out == NULL)
		error("Could not allocate memory\n");

	if(fftwf_init_threads())
		fftwf_plan_with_nthreads(CORES);

	plan1 = fftwf_plan_dft_2d(height, width, 
	                        img1, img1, FFTW_FORWARD, FFTW_ESTIMATE);

	plan2 = fftwf_plan_dft_2d(height, width, 
	                        img2, img2, FFTW_FORWARD, FFTW_ESTIMATE);

	planI = fftwf_plan_dft_2d(height, width, 
	                        img2, img2, FFTW_BACKWARD, FFTW_ESTIMATE);

	if(plan1 == NULL || plan2 == NULL || planI == NULL)
		error("Could not plan FFT\n");

	/* Initialize out array */
	for(int i = 0; i < px_count; i++) {
		out[i] = 1.0;
	}

	/* Calculate correlation for each band between the two images, multiply results of all bands */
	for(band = 1; band <= n_bands; band++) {
		printf("FFT 1 band %d\n", band);
		runFFT( plan1, srcDS1, img1, band, bbox1, width, height );
		printf("FFT 2 band %d\n", band);
		runFFT( plan2, srcDS2, img2, band, bbox2, width, height );

		
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
		//printf("Save band %d; min = %f max = %f\n", band, min, max);
		for(int i = 0; i < px_count; i++) {
			out[i] = out[i] * ((cabs(img2[i]) - min) / (max-min) );
		}

	}
	
	/* Cleanup memory */
	fftwf_destroy_plan(plan1);
	fftwf_destroy_plan(plan2);
	fftwf_destroy_plan(planI);
	fftwf_free(img1);
	fftwf_free(img2);

	fft2shift(out, width, height);

	/* Row and column of image center */
	const int pp_row = height / 2;
	const int pp_col = width / 2;

	const int height2 = pp_row * 2; /* Avoid odd pixels at end */
	const int width2 = pp_col * 2; /* image-cross-cor fft2shift has a bug with this */

	/* Find peak response that is outside the mask */
	float max = -1;
	int max_col = 0, max_row = 0;

	for(int row = 0; row < height2; row++) {
		int row_offset = row * width;

		for(int col = 0; col < width2; col++) {
			float v = out[row_offset + col];
			if(v > max) {
				max = v;
				max_col = col;
				max_row = row;
			}
		}
	}

	*offset_x = max_col - pp_col;
	*offset_y = max_row - pp_row;

	fftwf_free(out);

	return max;
}


void newBBox(int quadrant, const struct Rect* in, struct Rect* out)
{
	out->width = in->width / 2;
	out->height = in->height / 2;
	switch(quadrant) {
		case 1:
			out->minx = in->minx + (in->width)/2;
			out->miny = in->miny + (in->height)/2;
			break;
		case 2:
			out->minx = in->minx;
			out->miny = in->miny + (in->height)/2;
			break;
		case 3:
			out->minx = in->minx;
			out->miny = in->miny;
			break;
		case 4:
			out->minx = in->minx + (in->width)/2;
			out->miny = in->miny;
			break;
		default:
			out->minx = in->minx;
			out->miny = in->miny;
	}

}

int calcDisparityMap(GDALDataset* srcDS1, GDALDataset* srcDS2, GDALDataset* dstDS, const struct Rect bbox1, const struct Rect bbox2, const int initial_offset_x, const int initial_offset_y)
{
	struct Rect r1, r2;

	int offset_x = 0, offset_y = 0;
	int width, height, n_bands;

	int src_width, src_height;
	src_width = srcDS2->GetRasterXSize();
	src_height = srcDS2->GetRasterYSize();

	width = MY_MAX(bbox1.width, bbox2.width);
	height = MY_MAX(bbox1.height, bbox2.height);

	printf("Processing: (%d,%d) w: %d h: %d\n", bbox1.minx, bbox1.miny, bbox1.width, bbox1.height);

	float max = findOffset(srcDS1, srcDS2, bbox1, bbox2, width, height, &offset_x, &offset_y);

	offset_x += initial_offset_x;
	offset_y += initial_offset_y;

	printf("\t%d\t%d\t%f\n", offset_x, offset_y, max);

	if(bbox1.width > 64 && bbox1.height > 64 && bbox2.width > 64 && bbox2.height > 64 && max > 0.75) {
		// Quarter and try compute again.

		for(int i = 1; i <= 4; i++) {
			newBBox(i, &bbox1, &r1);	
			newBBox(i, &bbox2, &r2);

			/* Center search on expected area based on last try */
			r2.minx += offset_x; 
			r2.miny += offset_y;

			/* Make sure we don't go off the edge of the image */
			if(r2.minx < 0) {
			  offset_x = offset_x - r2.minx;
			  r2.minx = 0;
			}
			if(r2.miny < 0) {
			  offset_y = offset_y - r2.miny;
			  r2.miny = 0;
			}
			if(r2.minx + r2.width > src_width) {
			  offset_x = offset_x + r2.minx - (r2.width - src_width); //?
			  r2.minx =  r2.width - src_width;
			}
			if(r2.miny + r2.height > src_height) {
			  offset_y = offset_y + r2.miny - (r2.height - src_height); //?
			  r2.miny = r2.height - src_height;
			}

       			calcDisparityMap(srcDS1, srcDS2, dstDS, r1, r2, offset_x, offset_y);
		}
	} else {
		// Final Answer
		printf("(%d, %d) offset (%d, %d)\n", bbox1.minx + bbox1.width/2, bbox1.miny + bbox1.width/2, offset_x, offset_y);

		float disparity = sqrtf( offset_x * offset_x + offset_y * offset_y );
		dstDS->GetRasterBand(1)->RasterIO( GF_Write, bbox1.minx, bbox1.miny, bbox1.width, bbox1.height,
			   &disparity, 1, 1, GDT_CFloat32, 0, 0);

	}

	return 0;
}

int main(int argc, char** argv)
{
	char* srcFileName1;
	char* srcFileName2;
	char* dstFileName;

	GDALDataset* srcDS1, *srcDS2, *dstDS;
	
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

	dstDS = CreateOutputDataset(dstFileName, bbox1.width, bbox1.height, 1);

	calcDisparityMap(srcDS1, srcDS2, dstDS, bbox1, bbox2, 0, 0);

	delete srcDS1;
	delete srcDS2;
	delete dstDS;

	return 0;	
}
