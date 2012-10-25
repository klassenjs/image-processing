/*
g++ -I/Applications-POSIX/osgeo/include -I/opt/local/include -L/Applications-POSIX/osgeo/lib -L/opt/local/lib -lgdal -lproj -lfftw3f -lfftw3f_threads -lm -O2 image.cpp -o image

g++ -lgdal1.6.0 -lproj -lfftw3f -lfftw3f_threads -lm image3.cpp -o image3 -I/usr/include/gdal -L/home/jimk/PPP/local/lib -lopencv_core -lopencv_calib3d -lopencv_imgproc -lopencv_highgui -lopencv_contrib -I/home/jimk/PPP/local/include -Wl,-R -Wl,'/home/jimk/PPP/local/lib'

g++ -lgdal1.6.0 -lproj -lfftw3f -lfftw3f_threads -lm image3.cpp -o image3 -I/usr/include/gdal -L/home/jimk/PPP/local/lib -lopencv_core -lopencv_calib3d -lopencv_imgproc -lopencv_highgui -lopencv_contrib -I/home/jimk/PPP/local/include -Wl,-R -Wl,'/home/jimk/PPP/local/lib' -g
*/


/* Assumptions:
 *   Lens distortion removed
 *   Images rectified to each other
 *   Images overlap > 50%
 *   Most of image is valid (no large no-data areas)
 */

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <stdexcept>
#include <string>

/* For n-band Image IO */
#include <proj_api.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

/* For display and YAML */
#include <opencv2/core/core.hpp>
#include <opencv2/ts/ts.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>

/* For corelation = IFFT( FFT(A) * (i*FFT(B)) ) */ 
#include "/usr/include/complex.h"
#include <math.h>
#include <fftw3.h>

#include <getopt.h>

#ifndef CORES
 #define CORES (int)sysconf(_SC_NPROCESSORS_ONLN)
#endif

#define MY_MIN(a,b) (a<b)? a : b
#define MY_MAX(a,b) (a>b)? a : b

int gVERBOSE = 1;
int gWAITKEY = 10;

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

void info(const char* fmt, ... ) {
	va_list args;
	va_start(args, fmt);
	if(gVERBOSE) {
		vprintf(fmt, args);
	}
	va_end(args);
}


/***********************  calculateXYZ  ***********************/
float calculateXYZ2(
	cv::Mat ang1, cv::Mat t1p, cv::Mat ang2, cv::Mat t2p,
	float z1,
	float *X, float *Y, float *Z)
{
  cv::Mat world1 = cv::Mat::zeros(3,1,CV_32F);
  cv::Mat world2 = cv::Mat::zeros(3,1,CV_32F);

  world1 = (z1 * ang1) + t1p;

  // Solve for z2 to match Z1 and Z2 given z1    
  // world1 - t2 = z2 * ang2  |Z
  float z2 = (world1.at<float>(2) - t2p.at<float>(2)) / ang2.at<float>(2);

  world2 = (z2 * ang2) + t2p;

  float X1 = world1.at<float>(0);
  float Y1 = world1.at<float>(1);
  float X2 = world2.at<float>(0);
  float Y2 = world2.at<float>(1);
  float new_error = sqrt((X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2));

  //std::cout << "z1,err:" << z1 << "\t" << new_error << std::endl;    
  //std::cout << "world1: " << world1 << std::endl << "world2: " << world2 << std::endl << std::endl;

  *X = X1;
  *Y = Y1;
  *Z = world1.at<float>(2);

  return new_error;
}


/* calculateXYZ()
 * Inputs: 
 *    M1, M2 - Camera matricies
 *    R1, R2 - Rotation matrices
 *    t1, t2 - Translation vectors
 *    u1, u2 - column pixel coordinate
 *    v1, v2 - row pixel coordinate
 * Outputs:
 *    X, Y, Z - estimated position in 3D world space
 * Notes:
 *    x, y, z - position in camera coordinate system
 */
void calculateXYZ(const cv::Mat M1, const cv::Mat R1, const cv::Mat t1, 
		  const cv::Mat M2, const cv::Mat R2, const cv::Mat t2, 
		  int u1, int v1, int u2, int v2,
		  float* X, float* Y, float* Z)
{
  // Inverse the projection direction
  // OpenCV specifies [ camera_cs ] = R [ World_cs ] + t
  // We need to calculate world coordinates from camera coordinates
  // So Rnew = R'; tnew = R'(-t)

//  std::cout << "R1: " << R1 << std::endl << "R1': " << R1.inv() << std::endl;  
//  std::cout << "t1: " << t1 << std::endl << "t1': " << -R1.inv()*t1 << std::endl;  
//  std::cout << "R2: " << R2 << std::endl << "R2': " << R2.inv() << std::endl;  
//  std::cout << "t2: " << t2 << std::endl << "t2': " << -R2.inv()*t2 << std::endl;  

  // Need to copy results to new matrices otherwise args are modified.
  cv::Mat R1p = cv::Mat::zeros(3,3, CV_32F);
  cv::Mat R2p = cv::Mat::zeros(3,3, CV_32F);
  cv::Mat t1p = cv::Mat::zeros(3,1, CV_32F);
  cv::Mat t2p = cv::Mat::zeros(3,1, CV_32F);

  R1p = R1.inv();
  t1p = -R1p*t1;
  R2p = R2.inv();
  t2p = -R2p*t2;

  // Calculate the "angles" of the line extending from the focal center to the point on the ground
  float fx1 = M1.at<float>(0,0);  // Focal length
  float fy1 = M1.at<float>(1,1);
  float cx1 = M1.at<float>(0,2);  // Principle point
  float cy1 = M1.at<float>(1,2);

  float fx2 = M2.at<float>(0,0);
  float fy2 = M2.at<float>(1,1);
  float cx2 = M2.at<float>(0,2);
  float cy2 = M2.at<float>(1,2);

  cv::Mat ang1 = (cv::Mat_<float>(3,1) << 
	(u1 - cx1)/fx1,
	-(v1 - cy1)/fy1,
	-1 // fix for needed images swapped
  );
  cv::Mat ang2 = (cv::Mat_<float>(3,1) << 
	(u2 - cx2)/fx2,
	-(v2 - cy2)/fy2,
	-1 
  );  
  std::cout << "ang1: " << ang1 << std::endl << "ang2: " << ang2 << std::endl;  

  // Rotate into World CS
  ang1 = R1p * ang1;
  ang2 = R2p * ang2;

  std::cout << "R*ang1: " << ang1 << std::endl << "R*ang2: " << ang2 << std::endl;  


  // Note this is now a constrained optimization problem...  
  // ... minimize the difference between world1 and world2.
  // Overdetermined -- 3 eqns, 2 unknowns.  First reduce to 2 eqns, 1 unk.
  //  then LSQ.
  float X1, Y1, Z1;

  float z1 = 0.0; 
  float min_error = calculateXYZ2( ang1, t1p, ang2, t2p, z1, X, Y, Z );

  for(z1 = 1.0; z1 < t1p.at<float>(2); z1 += 4.0) {
    float error = calculateXYZ2( ang1, t1p, ang2, t2p, z1, &X1, &Y1, &Z1 );
    //std::cout << z1 << "\t" << error << std::endl;
    if(error < min_error) {
        min_error = error;
  	*X = X1;  
  	*Y = Y1;
  	*Z = Z1;
    }
  }
  return;
}

/* Globals - XYZ Calc */
bool skipXYZ = true;

cv::Mat M1, R1, t1;
cv::Mat M2, R2, t2;

void calculateXYZ(int u1, int v1, int u2, int v2,
		float *X, float* Y, float* Z)
{
  *X = 0.0; *Y = 0.0; *Z = 0.0;
  if(!skipXYZ) {
  	calculateXYZ(M1, R1, t1, M2, R2, t2, u1, v1, u2, v2, X, Y, Z);

  }
}


/***********************  calculate disparity  ***********************/

/*---------------------- FFT/Correlation Library ----------------------*/
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

/*--------------------------- Find offset of best match ---------------------*/

/* Globals - DisparityMap */
GDALDataset *srcDS1, *srcDS2, *dstDS;
FILE *xyzFile;
int   minSize = 32;
char* anaglyphBasename;
int   anaglyphMinSize = 512;

float findOffset(struct Rect bbox1, struct Rect bbox2, int width, int height, int* offset_x, int* offset_y) 
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
		info("FFT 1 band %d\n", band);
		runFFT( plan1, srcDS1, img1, band, bbox1, width, height );
		info("FFT 2 band %d\n", band);
		runFFT( plan2, srcDS2, img2, band, bbox2, width, height );

		
		info("Complex Conj band %d\n", band);
		/* mult img1 and conj of img2 */
		for(int px = 0; px < px_count; px++) {
			img2[px] = img1[px] * conj(img2[px]);
		}
	
		/* IFFT of result */	
		info("IFFT band %d\n", band);
		fftwf_execute(planI);	

		info("normalize band %d\n", band);
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
		//info("Save band %d; min = %f max = %f\n", band, min, max);
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


void displayMatch(struct Rect bbox1, struct Rect bbox2, int offset_x, int offset_y)
{
  /* Calculate size of output image */
  
  int minx = MY_MIN(bbox1.minx, bbox2.minx + offset_x);
  int miny = MY_MIN(bbox1.miny, bbox2.miny + offset_y);
  int maxx = MY_MAX(bbox1.minx+bbox1.width, bbox2.minx+bbox2.width+offset_x);
  int maxy = MY_MAX(bbox1.miny+bbox1.height, bbox2.miny+bbox2.height+offset_y);

  int width = maxx-minx;
  int height = maxy-miny;

  std::vector<cv::Mat> im;
  im.push_back( cv::Mat::zeros(height, width, CV_8U) );
  im.push_back( cv::Mat::zeros(height, width, CV_8U) );
  im.push_back( cv::Mat::zeros(height, width, CV_8U) );

  unsigned char* img1 = (unsigned char*)malloc(sizeof(unsigned char) * bbox1.width * bbox1.height);
  unsigned char* img2 = (unsigned char*)malloc(sizeof(unsigned char) * bbox1.width * bbox1.height);
  srcDS1->GetRasterBand(3)->RasterIO( GF_Read, bbox1.minx, bbox1.miny, 
				   bbox1.width, bbox1.height,
				   img1, bbox1.width, bbox1.height,
				   GDT_Byte, 0, 0);
  srcDS1->GetRasterBand(2)->RasterIO( GF_Read, bbox1.minx, bbox1.miny, 
				   bbox1.width, bbox1.height,
				   img2, bbox1.width, bbox1.height,
				   GDT_Byte, 0, 0);

  int dx = bbox1.minx - minx;
  int dy = bbox1.miny - miny;
  for(int y = 0; y < bbox1.height; y++) {
    for(int x = 0; x < bbox1.width; x++) {
      int px = x + y*bbox1.width;
      unsigned char px_val = img1[px];
      im[0].at<unsigned char>(y+dy, x+dx) = px_val;
      px_val = img2[px];
      im[1].at<unsigned char>(y+dy, x+dx) = px_val;
    }
  }
  free(img1);
  free(img2);

  img2 = (unsigned char*)malloc(sizeof(unsigned char) * bbox2.width * bbox2.height);
  srcDS2->GetRasterBand(1)->RasterIO( GF_Read, bbox2.minx, bbox2.miny, 
				   bbox2.width, bbox2.height,
				   img2, bbox2.width, bbox2.height,
				   GDT_Byte, 0, 0);

  dx = bbox2.minx+offset_x - minx;
  dy = bbox2.miny+offset_y - miny;
  for(int y = 0; y < bbox2.height; y++) {
    for(int x = 0; x < bbox2.width; x++) {
      int px = x + y*bbox2.width;
      unsigned char px_val = img2[px];
      im[2].at<unsigned char>(y+dy, x+dx) = px_val;
    }
  }
  free(img2);

  cv::Mat img = cv::Mat(height, width, CV_8UC3);
  cv::merge(im, img);
  cv::namedWindow("match", CV_WINDOW_NORMAL|CV_WINDOW_KEEPRATIO|CV_GUI_EXPANDED);
  cv::imshow("match", img);

  char saveName[1024];
  snprintf(saveName, 1024, "%s_%d-%d_%d.jpg", anaglyphBasename, bbox1.width, bbox1.minx, bbox1.miny);
  cv::imwrite(saveName, img);

  cv::waitKey(gWAITKEY);
}


// Learning search algorithm
int calcDisparityMap(const struct Rect bbox1, const struct Rect bbox2)
{
	struct Rect r1, r2;

	int offset_x = 0, offset_y = 0;
	int width, height, n_bands;

	int src_width, src_height;
	src_width = srcDS2->GetRasterXSize();
	src_height = srcDS2->GetRasterYSize();

	width = MY_MAX(bbox1.width, bbox2.width);
	height = MY_MAX(bbox1.height, bbox2.height);

	info("Processing: (%d,%d) w: %d h: %d\n", bbox1.minx, bbox1.miny, bbox1.width, bbox1.height);

	float max = findOffset(bbox1, bbox2, width, height, &offset_x, &offset_y);


	/* Find overlap between images */
	int initial_offset_x = bbox1.minx - bbox2.minx;
	int initial_offset_y = bbox1.miny - bbox2.miny;
	int total_offset_x = initial_offset_x + offset_x;
	int total_offset_y = initial_offset_y + offset_y;

	struct Rect overlap; /* in cs of srcDS1 */
	overlap.minx = MY_MAX(bbox1.minx, /*bbox2.minx +*/ total_offset_x);
	overlap.miny = MY_MAX(bbox1.miny, /*bbox2.miny +*/ total_offset_y);

	int maxx = MY_MIN(bbox1.minx + bbox1.width, total_offset_x + src_width /*bbox2.minx + total_offset_x + bbox2.width*/);
	int maxy = MY_MIN(bbox1.miny + bbox1.height, total_offset_y + src_height /*bbox2.miny + total_offset_y + bbox2.height*/);
	overlap.width = maxx - overlap.minx;
	overlap.height = maxy - overlap.miny;

	if(overlap.width <= 0 || overlap.height <= 0) {
	  error("no overlap");
	}

	struct Rect overlap2; /* in cs of srcDS2 */
       	overlap2.minx = overlap.minx - (total_offset_x);
	if(overlap2.minx < 0) {
	  overlap.minx = overlap.minx - overlap2.minx;
	  overlap.width = overlap.width - overlap2.minx;
	  overlap2.minx = 0;
	}

	overlap2.miny = overlap.miny - (total_offset_y);
	if(overlap2.miny < 0) {
	  overlap.miny = overlap.miny - overlap2.miny;
	  overlap.height = overlap.height - overlap2.miny;
	  overlap2.miny = 0;
	}

	maxx = srcDS2->GetRasterXSize() - (overlap2.minx + overlap.width);
	if(maxx < 0) {
	  overlap.width = overlap.width - maxx;
	}

	maxy = srcDS2->GetRasterYSize() - (overlap2.miny + overlap.height);
	if(maxy < 0) {
	  overlap.height = overlap.height - maxy;
	}
	
	overlap2.width = overlap.width;
	overlap2.height = overlap.height;

	info("overlap: %d %d %d %d\n", overlap.minx, overlap.miny, overlap.minx+overlap.width, overlap.miny+overlap.height);
	info("overlap2: %d %d %d %d\n", overlap2.minx, overlap2.miny, overlap2.minx+overlap2.width, overlap2.miny+overlap2.height);
	info("\t      offset: %d\t%d\t%f\n\ttotal offset: %d\t%d\n", offset_x, offset_y, max, total_offset_x, total_offset_y);

	/* Display Progress */
	if(overlap.width > anaglyphMinSize)
	  displayMatch(bbox1, bbox2, total_offset_x, total_offset_y);

	/* Termination condition */
	if(overlap.width > minSize && overlap.height > minSize && max > 0.75) {

		/* Operators -- Grid and compute again */
		int dim = MY_MIN(overlap.width, overlap.height);
		{ /* Force to be a power of two makes the following loop more efficient */
			unsigned char r = 0;
			while( dim >>= 1 )
				r++;
			dim = 1 << r;
		}	
		dim = (dim / 2) + (dim % 1);
		r1.width = r1.height = r2.width = r2.height = dim;

		for(int j = 0; j < overlap.height; j=j+dim) {
			int jj = j+dim > overlap.height ? overlap.height - dim : j;
			for(int i = 0; i < overlap.width; i=i+dim) {
				int ii = i+dim > overlap.width ? overlap.width - dim: i;
				r1.minx = overlap.minx + ii;
				r1.miny = overlap.miny + jj;
				r2.minx = overlap2.minx + ii;
				r2.miny = overlap2.miny + jj;
	       			calcDisparityMap(r1, r2);
			}
		}				
	} else {
	  /* Final Answer -- Compute XYZ */
          int u1 = overlap.minx + overlap.width/2;
          int v1 = overlap.miny + overlap.height/2;
          int u2 = overlap2.minx + overlap2.width/2;
          int v2 = overlap2.miny + overlap2.height/2;
	  info("(%d, %d) - (%d, %d)\n", u1, v1, u2, v2);

          /* Calculate XYZ */
	  if(!skipXYZ) {
          	float X,Y,Z;
          	calculateXYZ(u1, v1, u2, v2, &X, &Y, &Z);

	  	/* Save bands to XYZ file to colorize output */
          	int n_bands = srcDS1->GetRasterCount();
	  	char* bands;
          	int bands_len = 0;
	  	bands = (char*)alloca(10*n_bands);
	  	bzero(bands, 10*n_bands);
          	for(int band=1; band <= n_bands; band++) {
			float val;
		        srcDS1->GetRasterBand(band)->RasterIO( GF_Read, overlap.minx, overlap.miny,
        	                           overlap.width, overlap.height,
        	                           &val, 1, 1,
        	                           GDT_Float32, 0, sizeof(float));
			int chars = snprintf(bands + bands_len, 10*n_bands - bands_len, "%f ", val);
			if(chars > 0)
				bands_len += chars;
			else
				break;
		  }

	  	  fprintf(xyzFile, "%f %f %f %s\n", X, Y, Z, bands);
          }

          /* Save to disparity Map */
          if(dstDS) {
	  	float disparity = sqrtf( total_offset_x*total_offset_x + total_offset_y*total_offset_y );
	  	dstDS->GetRasterBand(1)->RasterIO( GF_Write, bbox1.minx, bbox1.miny, bbox1.width, bbox1.height,
						     &disparity, 1, 1, GDT_CFloat32, 0, 0);
          }
	}

	return 0;
}


GDALDataset* CreateOutputDataset(char* fileName, int width, int height, int bands)
{
	GDALDriverManager *gdm = GetGDALDriverManager();
	if(gdm == NULL)
		error("GetGDALDriverManager() failed!");

	GDALDriver *gd = gdm->GetDriverByName("GTiff");
	if(gd == NULL)
		error("Get GTiff Driver failed!");
	
	const char* ib = "INTERLEAVE=BAND"; // there must be a cleaner way...
	char ib2[20];
	strncpy(ib2, ib, 20);

	char* options[2];
	options[0] = ib2; //"INTERLEAVE=BAND";
	options[1] = NULL;
	
	GDALDataset *dstDS = gd->Create(fileName, width, height, bands, GDT_Float32, options);
	
	return dstDS;
}


void usage() {
	printf("image3 [options] left_image right_image\n");
	printf("\t--min-size=32              Minimum side length of pixel blocks to compare.\n");
	printf("\t--dispersion=[output_file] TIFF file to store dispersion map.\n");
	printf("\t--anaglyph=[basename]      Basename to store anaglyph snapshots. (e.g. basename_X_Y_size.png)\n");
	printf("\t--anaglyph-size=512        Smallest size anaglyph to save.\n");
	printf("\t--cores=[number]           Maximum number of threads to use during processing. (default #CPU cores)\n");
	printf("\n");
	printf("\t--xyz=[output_file]        ASCII XYZ file to store point cloud.\n");
	printf("\t--help                     Print this message and exit.\n");
}

int main(int argc, char** argv)
{
	char* srcFilename1 = NULL;
	char* srcFilename2 = NULL;
	char* dstFilename = NULL;
	char* xyzFilename = NULL;
	int cores = CORES;

	static const char *optString = "m:d:x:a:s:c:h";
	static const struct option longOpts[] = {
		{ "min-size", required_argument, NULL, 'm' },
		{ "dispersion", required_argument, NULL, 'd' },
		{ "xyz", required_argument, NULL, 'x' },
		{ "anaglyph", required_argument, NULL, 'a' },
		{ "anaglyph-size", required_argument, NULL, 's' },
		{ "cores", required_argument, NULL, 'c'},
		{ "help", no_argument, NULL, 'h'},
		{ NULL, no_argument, NULL, 0 }
	};
	int opt = 0;
	
	opt = getopt_long( argc, argv, optString, longOpts, NULL );
	while( opt != -1 ) {
		switch(opt) {
			case 'm':
				minSize = atoi( optarg );
				break;
			case 'd':
				dstFilename = optarg;
				break;
			case 'x':
				skipXYZ = false;
				xyzFilename = optarg;
				break;
			case 'a':
				anaglyphBasename = optarg;
				break;
			case 's':
				anaglyphMinSize = atoi( optarg );
				break;
			case 'c':
				cores = atoi( optarg );
				break;
			case 'h':
				usage();
				return(0);
				break;
			default:
				break;
		};
		opt = getopt_long( argc, argv, optString, longOpts, NULL );
	};

	if(argc != optind + 2) {
		usage();
		return(1);
	}
	srcFilename1 = argv[optind++];
	srcFilename2 = argv[optind++];

	info("Attempting to use %d cores for FFT\n", cores);
	if(fftwf_init_threads())
		fftwf_plan_with_nthreads(cores);

	GDALAllRegister();

	srcDS1 = (GDALDataset*) GDALOpen( srcFilename1, GA_ReadOnly );
	srcDS2 = (GDALDataset*) GDALOpen( srcFilename2, GA_ReadOnly );

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

	/* Initial search box is entire image */
	struct Rect bbox1, bbox2;
	bbox1.minx = bbox1.miny = bbox2.minx = bbox2.miny = 0;
	bbox1.width = srcDS1->GetRasterXSize();
	bbox1.height = srcDS1->GetRasterYSize();
	bbox2.width = srcDS2->GetRasterXSize();
	bbox2.height = srcDS2->GetRasterYSize();

	if(dstFilename)
		dstDS = CreateOutputDataset(dstFilename, bbox1.width, bbox1.height, 1);

	if(!skipXYZ) {
		xyzFile = fopen(xyzFilename, "w");

		cv::FileStorage fs1( (std::string(srcFilename1) + ".yaml").c_str(), cv::FileStorage::READ );
		fs1["A"] >> M1;
		fs1["R"] >> R1;
		fs1["T"] >> t1;
		fs1.release();
		cv::FileStorage fs2( (std::string(srcFilename2) + ".yaml").c_str(), cv::FileStorage::READ );
		fs2["A"] >> M2;
		fs2["R"] >> R2;
		fs2["T"] >> t2;
		fs2.release();

		M1.assignTo(M1, CV_32F);
		R1.assignTo(R1, CV_32F);
		t1.assignTo(t1, CV_32F);
		M2.assignTo(M2, CV_32F);
		R2.assignTo(R2, CV_32F);
		t2.assignTo(t2, CV_32F);
	}
	calcDisparityMap(bbox1, bbox2);

	delete srcDS1;
	delete srcDS2;
	delete dstDS;

	if(xyzFile)
		fclose(xyzFile);

	fftwf_cleanup_threads();

	return 0;	
}
