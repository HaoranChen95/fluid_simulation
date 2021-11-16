.PHONY: clean
VPATH=src include

ICPC = icpc -msse2 -qopt-report1 -qopt-report-phase=vec -qopt-report-phase=openmp -shared-intel -qopenmp -march=native -fp-model precise 
#GCC = g++ -msse2 -ffast-math -Wall -g -march=native 
GCC = g++ -msse2 -fopenmp -ffast-math -Wall -g -march=native
# GCC = g++ -msse2 -fopenmp -ffast-math -Wall -g -march=native -Q --help=target -v
CFLAGS = -std=c++0x -O2 -mcmodel=medium
# CFLAGS = -std=c++17 -O2 -mcmodel=medium
objs = main.o

simple_fluid_simulation: initialization.o implement.o export.o main.o
	$(ICPC) $(CFLAGS) -o $@ $^ -I include $(LDFLAGS)

%.o: %.cpp %.hpp
	$(ICPC) $(CFLAGS) -c -o $@ $^ -I include $(LDFLAGS)

clean:
	rm -rf *.o *.optrpt sf_simu
