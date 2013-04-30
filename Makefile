all: pyfinterpol.so pyfsla.so
pyfinterpol.so: pyfinterpol.f90
	f2py2 -c -m pyfinterpol{,.f90}
pyfsla.so: pyfsla.f90
	f2py2 -c -m pyfsla{,.f90} -lsla
pyfwcs.so: pyfwcs.f90
	f2py2 -I. -c -m pyfwcs{,.f90} -lwcs

clean:
	rm -rf *.so *.pyc *.o __pycache__
