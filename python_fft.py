import Image, numpy

def calcOffset(img1, img2, imgO):
	print "Loading " + img1
	i = Image.open(img1)
	a = numpy.asarray(i)
	f = numpy.fft.fft2(a)
	a = None

	print "Loading " + img2
	i2 = Image.open(img2)
	a2 = numpy.asarray(i2)
	f2 = numpy.fft.fft2(a2)
	a2 = None

	print "Calculating cross-correlation"
	f2c = f2.conj()
	f2 = None

	fi = f * f2c
	f = None
	f2c = None

	ai = numpy.fft.ifft2(fi)
	fi = None

	print "Calculating Offset"
	ar = numpy.real(ai)
	
	# Normalize Array (0 - 1)
	ar = ar - numpy.min(ar)
	ar = ar * (1/numpy.max(ar))
	
	ar_max = numpy.max(ar)
	print numpy.nonzero(ar > (.7 * ar_max))


	print "Saving output results as " + imgO
	ir = Image.fromarray(ar)
	ir.save(imgO)

	# Note missing FFT shift...

calcOffset("56038_158_1.tif", "56038_159_1.tif", "out.tif")
# Max at 3172, 45 -> 159 is 3179px east, 45px south of 158
calcOffset("56038_159_1.tif", "56038_160_1.tif", "out2.tif")
# Max at 3200, 169 -> 6372px east, 169px south of 158
calcOffset("56038_158_1.tif", "56038_160_1.tif", "out3.tif")
# Max is invalid (there is a local max near where it should be,
# but it looks like offsets > 50% of the image cause issues)
