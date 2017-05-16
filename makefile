# call : make "name of program"
FC=gfortran
objects=beamform_ncf.o clogc.o sacio.o
executable=beamforming
all: sacio.mod $(executable) 
%.o:%.f90
	$(FC) $(FFLAG1) $^ -c
sacio.mod:sacio.f90
	$(FC) $^ -c
$(executable):$(objects)
	$(FC) $^ -o $@
install:
	cp $(executable) ../bin/
uninstall:
	-rm ~/bin/$(excecutable)
clean:
	-rm $(objects) 
