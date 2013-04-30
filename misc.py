import numpy as np, os, errno, h5py as h5
import pyfftw

class Bunch(object):
	def __init__(self, **kwds):
		self.__dict__.update(kwds)
	def __eq__(self, other):
		if type(self) != type(other): return False
		return self.__dict__ == other.__dict__

def ctime2mjd(ctime):
	return ctime/86400 + 40587.0
deg2rad = np.pi/180
am2rad = deg2rad/60
as2rad = am2rad/60
degree = deg2rad
arcmin = am2rad
arcsec = as2rad

def mkdir_safe(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise

# Given a and b, both logically
# cyclic with period p, return
# recentered version of b such
# that each point in b is at most
# +- p/2 away from the corresponding
# point in a
def recycle(a,b,p): return a + (b-a+p/2) % p - p/2

#myfft = np.fft
myfft = pyfftw.interfaces.numpy_fft
pyfftw.interfaces.cache.enable()

# To keep fourier space units independent of the length
# of the array, we will work with (1/n**0.5,1/n**0.5)-normalized
# ffts. This differs from numpy, which uses (1,1/n) as the
# normalization. The latter is good for convolutions and
# decimation, but not for noise models.
def rfft(a):
	n = a.shape[-1]
	b = a.reshape(np.prod(a.shape[:-1]),n)
	fb = np.empty((b.shape[0],b.shape[1]/2+1),dtype=np.result_type(a,0j))
	for i in range(b.shape[0]):
		fb[i] = myfft.rfft(b[i])*n**-0.5
	return np.reshape(fb, list(a.shape[:-1]) + [fb.shape[-1]])
# Numpy irfft divides by 1/n, but we have already divided
# by 1/n**0.5, so we need to compensate here
def irfft(fa,n):
	fb = fa.reshape(np.prod(fa.shape[:-1]),fa.shape[-1])
	b  = np.empty((fb.shape[0],n))
	for i in range(b.shape[0]):
		b[i] = myfft.irfft(fb[i],n)*n**0.5
	return np.reshape(b, list(fa.shape[:-1]) + [b.shape[-1]])

def rfft_inplace(a, b):
	for d in range(a.shape[0]):
		b[d] = myfft.rfft(a[d])*a.shape[-1]**-0.5

def irfft_inplace(a, b):
	for d in range(a.shape[0]):
		b[d] = myfft.irfft(a[d])*b.shape[-1]**0.5

def ffunion(shape, dtype=np.float64):
	buf = pyfftw.n_byte_align_empty(list(shape[:-1])+[(shape[-1]/2+1)*2],16,dtype=dtype)
	#buf = np.empty(list(shape[:-1])+[(shape[-1]/2+1)*2],dtype=dtype)
	tod = buf[...,:shape[-1]]
	ft  = buf.view(dtype=np.result_type(dtype,0j))
	return tod, ft

def parse_slice(desc):
	class Foo:
		def __getitem__(self, p): return p
	foo = Foo()
	return eval("foo"+desc)

def h5dump(fname, data):
	mkdir_safe(os.path.dirname(fname))
	with h5.File(fname,"w") as hfile:
		hfile["data"] = data

def h5get(fname, field):
	with h5.File(fname,"r") as hfile:
		return hfile[field]

reset   = "\033[0m"
black   = "\033[0;30m"
red     = "\033[0;31m"
green   = "\033[0;32m"
brown   = "\033[0;33m"
blue    = "\033[0;34m"
purple  = "\033[0;35m"
cyan    = "\033[0;36m"
lgray   = "\033[0;37m"
gray    = "\033[1;30m"
lred    = "\033[1;31m"
lgreen  = "\033[1;32m"
lbrown  = "\033[1;33m"
lblue   = "\033[1;34m"
lpurple = "\033[1;35m"
lcyan   = "\033[1;36m"
white   = "\033[1;37m"
