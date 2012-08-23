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




void displayMatch(GDALDataset* srcDS1, GDALDataset* srcDS2, struct Rect bbox1, struct Rect bbox2, int offset_x, int offset_y)
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
  snprintf(saveName, 1024, "snapshot_%d-%d_%d.jpg", bbox1.width, bbox1.minx, bbox1.miny);
  cv::imwrite(saveName, img);

  cv::waitKey(gWAITKEY);
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
  float cx1 = M1.at<float>(2,0);  // Principle point
  float cy1 = M1.at<float>(2,1);

  float fx2 = M2.at<float>(0,0);
  float fy2 = M2.at<float>(1,1);
  float cx2 = M2.at<float>(2,0);
  float cy2 = M2.at<float>(2,1);

  cv::Mat ang1 = (cv::Mat_<float>(3,1) << 
	(u1 - cx1)/fx1,
	(v1 - cy1)/fy1,
	1 
  );
  cv::Mat ang2 = (cv::Mat_<float>(3,1) << 
	(u2 - cx2)/fx2,
	(v2 - cy2)/fy2,
	1 
  );  
//  std::cout << "ang1: " << ang1 << std::endl << "ang2: " << ang2 << std::endl;  

  // Rotate into World CS
  ang1 = R1p * ang1;
  ang2 = R2p * ang2;

//  std::cout << "R*ang1: " << ang1 << std::endl << "R*ang2: " << ang2 << std::endl;  

  cv::Mat world1 = cv::Mat::zeros(3,1,CV_32F);
  cv::Mat world2 = cv::Mat::zeros(3,1,CV_32F);
  float z1, z2;

  // Note this is now a constrained optimization problem...  
  // ... minimize the difference between world1 and world2.
  // Overdetermined -- 3 eqns, 2 unknowns.  First reduce to 2 eqns, 1 unk.
  //  then LSQ.

  float error = -1;
  for(z1 = 0.0; z1 > -10000.0; z1 -= 1.0) {
    world1 = (z1 * ang1) + t1p;

    // Solve for z2 to match Z1 and Z2 given z1    
    // world1 - t2 = z2 * ang2  |Z
    z2 = (world1.at<float>(2) - t2p.at<float>(2)) / ang2.at<float>(2);

    world2 = (z2 * ang2) + t2p;

//    std::cout << "z1: " << z1 << std::endl;
//    std::cout << "world1: " << world1 << std::endl << "world2: " << world2 << std::endl << std::endl;

    float X1 = world1.at<float>(0);
    float Y1 = world1.at<float>(1);
    float X2 = world2.at<float>(0);
    float Y2 = world2.at<float>(1);
    float new_error = (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2);
//    std::cout << "err: " << new_error << std::endl;    

    // Best so far...
    if(error < 0 || new_error < error) {
        *X = X1;
        *Y = Y1;
        *Z = world1.at<float>(2);
	error = new_error;
    }
  }
  return;
}

void testXYZ2(int u1, int v1, int u2, int v2)
{
  cv::Mat M1 = (cv::Mat_<float>(3,3) << 10000.0, 0.0, 3839.5, 0.0, 10000.0, 6911.5, 0.0, 0.0, 1.0);
  cv::Mat R1 = (cv::Mat_<float>(3,3) << -0.73056234718734, -0.6828392198637, 0.003042480576291, 0.68282529346054, -0.73049824782602, 0.011042125413153, -0.0053174695727292, 0.010144443752124, 0.99993440523782);
  cv::Mat t1 = (cv::Mat_<float>(3,1) << 3946076.5220691, 3398275.7654668, -52833.142739149);

  cv::Mat M2 = (cv::Mat_<float>(3,3) << 10000.0, 0.0, 3839.5, 0.0, 10000.0, 6911.5, 0.0, 0.0, 1.0);
  cv::Mat R2 = (cv::Mat_<float>(3,3) << -0.73024011297182, -0.68315173576369, 0.0072858307670647, 0.68318691406176, -0.7301596949777, 0.011066177516865, -0.0022400584083765, 0.013058550958209, 0.99991222434032);
  cv::Mat t2 = (cv::Mat_<float>(3,1) << 3948436.9257857, 3396314.5768609, -69657.736900173);

  float X,Y,Z;

  calculateXYZ(M1, R1, t1, M2, R2, t2, u1, v1, u2, v2, &X, &Y, &Z);

  printf("FOUND: %f %f %f\n", X, Y, Z);
}


void testXYZ()
{
  cv::Mat M1 = (cv::Mat_<float>(3,3) << 10000.0, 0.0, 3839.5, 0.0, 10000.0, 6911.5, 0.0, 0.0, 1.0);
  cv::Mat R1 = (cv::Mat_<float>(3,3) << -0.73056234718734, -0.6828392198637, 0.003042480576291, 0.68282529346054, -0.73049824782602, 0.011042125413153, -0.0053174695727292, 0.010144443752124, 0.99993440523782);
  cv::Mat t1 = (cv::Mat_<float>(3,1) << 3946076.5220691, 3398275.7654668, -52833.142739149);

  cv::Mat M2 = (cv::Mat_<float>(3,3) << 10000.0, 0.0, 3839.5, 0.0, 10000.0, 6911.5, 0.0, 0.0, 1.0);
  cv::Mat R2 = (cv::Mat_<float>(3,3) << -0.73024011297182, -0.68315173576369, 0.0072858307670647, 0.68318691406176, -0.7301596949777, 0.011066177516865, -0.0022400584083765, 0.013058550958209, 0.99991222434032);
  cv::Mat t2 = (cv::Mat_<float>(3,1) << 3948436.9257857, 3396314.5768609, -69657.736900173);

  float X,Y,Z;
  int u1,v1,u2,v2;
  u1 = 0;
  v1 = 0;
  u2 = -3000;
  v2 = 0;

  std::cout << "R1 was: " << R1 << std::endl;
  calculateXYZ(M1, R1, t1, M2, R2, t2, u1, v1, u2, v2, &X, &Y, &Z);
  std::cout << "R1 is: " << R1 << std::endl;
  printf("FOUND: %f %f %f\n", X, Y, Z);
}

void newBBox(int quadrant, const struct Rect* in, struct Rect* out)
{
	out->width = in->width / 2;
	out->height = in->height / 2;

	/* Avoid Gaps */
	if(2 * out->width < in->width)
		out->width++;
	if(2 * out->height < in->height)
		out->height++;

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

// Learning search algorithm
int calcDisparityMap(GDALDataset* srcDS1, GDALDataset* srcDS2, GDALDataset* dstDS, const struct Rect bbox1, const struct Rect bbox2)
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

	float max = findOffset(srcDS1, srcDS2, bbox1, bbox2, width, height, &offset_x, &offset_y);


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

	if(overlap.width > 512)
	  displayMatch(srcDS1, srcDS2, bbox1, bbox2, total_offset_x, total_offset_y);

	if(overlap.width > 32 && overlap.height > 32 && max > 0.75) {

		// Grid and compute again
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
	       			calcDisparityMap(srcDS1, srcDS2, dstDS, r1, r2);
			}
		}				
/*
	  	// Quarter and try compute again.
		for(int i = 1; i <= 4; i++) {
			newBBox(i, &overlap, &r1);	
			newBBox(i, &overlap2, &r2);

       			calcDisparityMap(srcDS1, srcDS2, dstDS, r1, r2);
		}
*/
	} else {
	  // Final Answer
	  int dx = overlap.minx + overlap.width/2 - overlap2.minx + overlap2.width/2;
	  int dy = overlap.miny + overlap.height/2 - overlap2.miny + overlap2.height/2;
	  printf("(%d, %d) -> (%d, %d)\n", overlap.minx + overlap.width/2, overlap.miny + overlap.height/2,
		  overlap2.minx + overlap2.width/2, overlap2.miny + overlap2.height/2);

          testXYZ2(overlap.minx + overlap.width/2, overlap.miny + overlap.height/2, 
                   overlap2.minx + overlap2.width/2, overlap2.miny + overlap2.height/2);

	  float disparity = sqrtf( total_offset_x*total_offset_x + total_offset_y*total_offset_y );
	  dstDS->GetRasterBand(1)->RasterIO( GF_Write, bbox1.minx, bbox1.miny, bbox1.width, bbox1.height,
					     &disparity, 1, 1, GDT_CFloat32, 0, 0);

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

int main(int argc, char** argv)
{
	char* srcFileName1;
	char* srcFileName2;
	char* dstFileName;

	//testXYZ();
	//return(0);

	info("Attempting to use %d cores for FFT\n", CORES);
	if(fftwf_init_threads())
		fftwf_plan_with_nthreads(CORES);

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

	calcDisparityMap(srcDS1, srcDS2, dstDS, bbox1, bbox2);

	delete srcDS1;
	delete srcDS2;
	delete dstDS;

	fftwf_cleanup_threads();

	return 0;	
}
