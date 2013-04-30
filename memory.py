# Get information about the process's memory usage.
import resource, os

def current():
	with open("/proc/%d/statm" % os.getpid(),"r") as f:
		return int(f.readline().split()[0])*resource.getpagesize()

def max():
	with open("/proc/%d/status" % os.getpid(),"r") as f:
		for line in f:
			toks = line.split()
			if toks[0] == "VmPeak:":
				return int(toks[1])*1024

def dmem(str):
	print "%10.4f %10.4f %s" % (current()/1024.**3, max()/1024.**3, str)
