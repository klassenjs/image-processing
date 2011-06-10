
# MacOS X
#CPPFLAGS=-I/Applications-POSIX/osgeo/include -I/opt/local/include -DCORES=2
#LDFLAGS=-L/Applications-POSIX/osgeo/lib -L/opt/local/lib -lgdal -lproj -lfftw3f -lfftw3f_threads -lm 

# Linux
CPPFLAGS=-I/usr/include/gdal -DCORES=8 -O3 -march=native
LDFLAGS=-lgdal1.6.0 -lproj -lfftw3f -lfftw3f_threads -lm

all: image-cross-cor image-fft complex2rgb image-find-cpp

image-cross-cor: image-cross-cor.cpp

image-fft: image-fft.cpp

complex2rgb: complex2rgb.cpp

image-find-cpp: image-find-cpp.cpp
