/*
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

void error(const char* msg)
{
	throw(std::out_of_range(msg));
} 


int CalcFFT(GDALDataset *srcDS, GDALDataset *dstDS) 
{
	fftwf_plan plan;
	fftwf_complex *in;  /* This is a complex double (64-bit real + 64-bit imag) */
	fftwf_complex *out;
	int band;

	/* TODO: Does GDAL give us a copy or does it own it... if we get a copy
  	 *       then we can work in place (out = in)
	 */
	const size_t px_count = srcDS->GetRasterXSize() * srcDS->GetRasterYSize();
	const size_t buffer_len = sizeof(fftwf_complex) * px_count;
	
	printf("Allocating %ld MB for FFT\n", buffer_len / (1024 * 1024));
	in  = (fftwf_complex*) fftwf_malloc(buffer_len);
	out = in; //(fftwf_complex*) fftw_malloc(buffer_len);

	if(fftwf_init_threads())
		fftwf_plan_with_nthreads(2);

	plan = fftwf_plan_dft_2d(srcDS->GetRasterYSize(), srcDS->GetRasterXSize(), 
	                        in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for(band = 1; band <= srcDS->GetRasterCount(); band++) {	
		printf("Processing band %d\n", band);
		srcDS->GetRasterBand(band)->RasterIO( GF_Read, 0, 0, 
	                                   srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
	                                   in, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
	                                   GDT_CFloat32, 0, 0);
	
		for(int i = 0; i < 10; i++) {
			printf("%f + %fj ", creal(in[i]), cimag(in[i]));
		}
		printf("\n");
		
		printf("Running FFT\n");
		fftwf_execute(plan);

		printf("Normalizing FFT\n");
		complex double norm = csqrt(px_count + 0I);
		for(int i = 0; i < px_count; i++) {
			out[i] = out[i] / norm;
		}
		
		for(int i = 0; i < 10; i++) {
			printf("%f + %fj ", creal(out[i]), cimag(out[i]));
		}
		printf("\n");
		
		printf("Saving Result\n");
		dstDS->GetRasterBand(band)->RasterIO( GF_Write, 0, 0,
	                                   srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
	                                   out, srcDS->GetRasterXSize(), srcDS->GetRasterYSize(),
	                                   GDT_CFloat32, 0, 0);
		
	}

	fftwf_destroy_plan(plan);
	fftwf_free(in);
	//fftw_free(out);
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
	
	GDALDataset *dstDS = gd->Create(fileName, width, height, bands, GDT_CFloat32, options);
	
	return dstDS;
}


int main(int argc, char** argv)
{
	char* srcFileName;
	char* dstFileName;

	GDALDataset* srcDS;
	GDALDataset* dstDS;

	if(argc != 3)
		error("Usage: prog infile outfile");

	srcFileName = argv[1];
	dstFileName = argv[2];

	GDALAllRegister();

	srcDS = (GDALDataset*) GDALOpen( srcFileName, GA_ReadOnly );

	if(srcDS == NULL)
		error("Could not open source dataset");

	printf("Got image %d x %d x %d\n",
	      srcDS->GetRasterXSize(),
	      srcDS->GetRasterYSize(),
	      srcDS->GetRasterCount()
	);

	dstDS = CreateOutputDataset(dstFileName, 
	                            srcDS->GetRasterXSize(),
	                            srcDS->GetRasterYSize(),
	                            srcDS->GetRasterCount()
	                           );

	if(dstDS == NULL)
		error("Could not create destination dataset");

	CalcFFT(srcDS, dstDS);

	delete srcDS;
	delete dstDS;

	return 0;	
}
