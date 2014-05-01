# Copyright (c) 2011-2013, James Klassen
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# MacOS X
#CPPFLAGS=-I/Applications-POSIX/osgeo/include -I/opt/local/include -DCORES=2
#LDFLAGS=-L/Applications-POSIX/osgeo/lib -L/opt/local/lib -lgdal -lproj -lfftw3f -lfftw3f_threads -lm 

# Linux
CXX=clang++

CPPFLAGS=-I/usr/include/gdal -DCORES=8
CXXFLAGS=-O4 -march=native -mtune=native -g
LDLIBS=-lgdal1.7.0 -lproj -lfftw3f -lfftw3f_threads -lm

all: image-cross-cor image-fft complex2rgb image-find-cpp image2 image3

image3: image3.cpp
	$(CXX) $(CXXFLAGS) \
		image3.cpp -o image3 \
		-I/usr/include/gdal \
		-I/home/jimk/apps/OpenCV/2.4.2/include \
		-L/home/jimk/apps/OpenCV/2.4.2/lib \
		-Wl,-R -Wl,'/home/jimk/apps/OpenCV/2.4.2/lib' \
		$(LDLIBS) \
		-lopencv_core -lopencv_calib3d -lopencv_imgproc \
		-lopencv_highgui -lopencv_contrib

