FC         = gfortran

FFLAGS = -O3 -g -pedantic -C

AR  = ar
RM  = -rm -v
SRC     := $(wildcard *.f)
OBJECTS := $(SRC:.f=.o)

all : ncosmo3.a test_ncosmo3 

ncosmo3.a : $(OBJECTS)
	$(FC) -c -g -pedantic -O3 $(SRC)
	$(AR) -rv ncosmo3.a $(OBJECTS)

test_ncosmo3 : test_ncosmo3.o
	$(FC) $(FFLAGS) -o test_ncosmo3 \
	ncosmo3.a

.PHONY: clean
clean:
	$(RM) $(OBJECTS) ncosmo3.a test_ncosmo3
