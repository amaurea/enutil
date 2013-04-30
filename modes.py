# Solve for the common mode and the individual amplitudes.
# Data model is d[d,s] = A[m,d]*common[m,s] + rest[d,s].
# m is the number of common modes to look for. For temperature
# this would normally be 1, and for polarizaiton it would be
# 3.
import numpy as np
from scipy.optimize import fmin_ncg, fmin_powell

# This solution does not work for m > 1 due to a degeneracy between
# the various modes.
def find_common_mode(d, m):
	# Set up problem set. Assume 1 as starting value for parameters
	info = Info(d, m)
	# Set up a reasonable starting guess: All amplitudes are 1. Then
	# do a simple amps->mode->amps iteration to handle the few amplitudes
	# that are reversed.
	amps  = np.zeros((info.nmode, info.ncol))+1
	mode  = mode_from_amps(d, amps)
	amps  = amps_from_mode(d, mode)

	# Solve the system
	param = flatten(amps)
	#param = fmin_ncg(chisq, param, dchisq, args=(info,), disp=False, maxiter=16)
	param = fmin_powell(chisq, param, args=(info,), disp=False, maxiter=4)
	#param = minimize(chisq, param, args=(info,)).x
	# Extract the solution
	amps = expand(param, info.nmode, info.ncol)
	mode = mode_from_amps(d, amps)
	return mode, amps

class Info:
	def __init__(self, d, m):
		self.d = d
		self.nsamp = d.shape[1]
		self.ncol  = d.shape[0]
		self.nmode = m
		self.norm  = 1.0/np.sum(d**2)
		self.mode = None
		self.amps = None
	def get_mode(self, amps):
		if self.amps == None or np.any(amps != self.amps):
			self.amps = np.array(amps)
			self.mode = mode_from_amps(self.d, amps)
		return self.mode

def chisq(param, info):
	amps = expand(param, info.nmode, info.ncol)
	mode = info.get_mode(amps)
	residual = info.d-amps.T.dot(mode)
	res = 0.5*np.trace(residual.dot(residual.T))*info.norm
	#print "In chisq", np.min(param), np.max(param), np.mean(param), res
	return res

def dchisq(param, info):
	amps = expand(param, info.nmode, info.ncol)
	mode = info.get_mode(amps)
	residual = info.d-amps.T.dot(mode)
	res = flatten(-mode.dot(residual.T))*info.norm
	return res

def dchisq_emp(param, info):
	res = np.empty((len(param)))
	delta = 1e-5
	for d in range(len(res)):
		p1, p2 = np.array(param), np.array(param)
		p1[d] += delta
		p2[d] -= delta
		res[d] = (chisq(p1, info)-chisq(p2,info))/(2*delta)
	return res

def mode_from_amps(d, amps):
	# mode = d dot amps dot (amps.T dot amps)**-1
	# but indices are transposed here compared to my calculation
	A = amps.dot(amps.T)
	b = amps.dot(d)
	return np.linalg.solve(A,b)

def amps_from_mode(d, mode):
	A = mode.dot(mode.T)
	b = mode.dot(d.T)
	return np.linalg.solve(A,b)

def flatten(amps):
	# Go from ndet by m amps to (ndet-1)*m flattened amps.
	# The -1 is due to the degeneracy between mode scaling and
	# the amps.
	tmp = amps[:,1:]
	return np.reshape(tmp, (tmp.size))

def expand(param, nmode, ncol):
	tmp = np.empty((nmode,ncol))
	tmp[:,0] = 1
	tmp[:,1:] = param.reshape(nmode,ncol-1)
	return tmp
