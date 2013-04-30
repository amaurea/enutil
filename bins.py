# Various binnin schemes
import numpy as np

def linbin(n, nbin=None):
	if not nbin: nbin = int(np.round(n**0.5))
	tmp  = np.arange(nbin+1)*n/nbin
	return np.vstack((tmp[:-1],tmp[1:])).T

def expbin(n, nbin=None, nmin=8):
	if not nbin: nbin = int(np.round(n**0.5))
	tmp  = np.array(np.exp(np.arange(nbin+1)*np.log(n+1)/nbin)-1,dtype=int)
	fixed = [tmp[0]]
	i = 0
	while i < nbin:
		for j in range(i+1,nbin+1):
			if tmp[j]-tmp[i] >= nmin:
				fixed.append(tmp[j])
				i = j
	tmp = np.array(fixed)
	tmp[-1] = n
	return np.vstack((tmp[:-1],tmp[1:])).T
