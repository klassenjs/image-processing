
# MacOS X
#CPPFLAGS=-I/Applications-POSIX/osgeo/include -I/opt/local/include -DCORES=2
#LDFLAGS=-L/Applications-POSIX/osgeo/lib -L/opt/local/lib -lgdal -lproj -lfftw3f -lfftw3f_threads -lm 

# Linux
CPPFLAGS=-I/usr/include/gdal -DCORES=8
CXXFLAGS=-Os -march=native -mtune=native
LDLIBS=-lgdal1.7.0 -lproj -lfftw3f -lfftw3f_threads -lm

all: image-cross-cor image-fft complex2rgb image-find-cpp image2 image3

image3: image3.cpp
	g++ \
                $(CXXFLAGS) \
		image3.cpp -o image3 \
		-I/usr/include/gdal \
		-I/home/jimk/apps/OpenCV/2.4.2/include \
		-L/home/jimk/apps/OpenCV/2.4.2/lib \
		-Wl,-R -Wl,'/home/jimk/apps/OpenCV/2.4.2/lib' \
		-lgdal1.7.0 -lproj \
		-lfftw3f -lfftw3f_threads -lm \
		-lopencv_core -lopencv_calib3d -lopencv_imgproc \
		-lopencv_highgui -lopencv_contrib

