# This module provides pointing interpolation, to allow rapid calculation of
# pointing for all detectors on the fly. This is needed because the full
# pointing is too large to keep in memory during the CG chain.
#
# In my use, the interpolator will turn ([mjd,-az,pi/2-el],samp) into
# ([ra,pi/2-dec,cos2psi,sin2psi],samp). It is a simple 3d grid in the input,
# and returns a ax+b, where x is the distance from the closest grid point,
# and a and b are stored in that point.
#
# The interpolator will store the following quantities:
#  self.y: (ngrid,odim), self.dy: (ngrid,odim,idim), self.box(2,idim),
#  self.n: (idim)
# The inputs and outputs will have the ordering:
#  x: (nsamp,idim), y: (nsamp,odim)
# y and dy are logically N-dimensional grids, but variable dimensionality is
# not practical to work with.
import numpy as np, pyfinterpol

class NdInterpol:
	def __init__(self, box, n, vals):
		self.box  = np.array(box,dtype=float)
		self.n    = np.array(n)
		self.y    = np.array(vals)
		self.dy   = calc_gradients(self.y, self.n+1)
	# Return interpolated values at every position in x.
	def __call__(self, x):
		return pyfinterpol.ipol(x.T, self.box.T, self.n+1, self.y.T, self.dy.T).T

# Builds an interpolation to required precision.
# This implementation is somewhat wasteful, in
# that it evaluates some items more often than required.
# The alternative would, however, be significantly less
# clear.
def build_interpol(box, func, maxerr):
	box        = np.array(box,dtype=float)
	idim       = box.shape[1]
	n          = np.array([3]*idim)
	xgrid      = makegrid(box[0], box[1], n)

	# Set up initial interpolation
	ip = NdInterpol(box, n, func(xgrid))

	# Refine until good enough
	errs = [None]*idim
	while True:
		nok = 0
		# Consider accuracy for each input parameter
		for i in range(idim):
			if errs[i] == None or any(errs[i] > maxerr):
				# Grid may not be good enough in this direction.
				# Try doubling resolution
				nnew = np.array(ip.n); nnew[i] *= 2
				xnew = makegrid(box[0], box[1], nnew)
				yinter = ip(xnew)
				ytrue  = func(xnew)
				diff = ytrue-yinter
				ind  = np.argmax(diff,0)
				err = np.amax(abs(ytrue-yinter), 0)
				if any(err > maxerr):
					# Not good enough, so accept improvement
					ip = NdInterpol(box, nnew, ytrue)
				else: nok += 1
				errs[i] = err
			else: nok += 1
		if nok >= idim: break
	return ip

# Return an (ngrid(n+1),idim) array of evenly spaced points
# within the rectangle specified by x0:x1. The mapping from
# the implicit N-dimensional space to the actually used 1-d
# space is the same as the one used in the interpolation method.
# n is the number of cells - the number of points in each direction
# is n+1.
def makegrid(x0, x1, n):
	inds  = np.rollaxis(np.indices(n+1),0,len(n)+1) # (d1,d2,d3,...,indim)
	finds = np.reshape(inds, (np.prod(inds.shape[:-1]), inds.shape[-1]))
	return finds * (x1-x0)/n + x0

# Given y(ngrid,odim) and n(idim) return the gradient dy(ngrid,odim,idim).
# Implemented by expanding to a full nD array, calling numpy's gradient,
# and then collapsing back.
def calc_gradients(y, n):
	odim  = y.shape[-1]
	idim  = len(n)
	nd_y  = np.reshape(y, list(n) + [odim])
	grads = np.array([np.reshape(np.gradient(nd_y[...,i]),(idim,np.prod(n))) for i in range(odim)])
	# Result is now (odim,idim,ngrid), so move ngrid to 0
	return np.ascontiguousarray(np.rollaxis(grads,2))
