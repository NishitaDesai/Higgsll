CERNLIB = -Wl,-static -lpdflib804 -lmathlib \
      -lpacklib -lkernlib -Wl,-dy	\
      -lm -lnsl -lcrypt -ldl -lgfortran
LIBLHA = -L /data/lhapdf/lib -lLHAPDF
LIB = -lgfortran -lgcc	
DEBUG = -O0 -ggdb3
#CERNPATHS = -static -L/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.4/x86_64-slc5-gcc43-opt/lib -lLHAPDF \
	-lgfortran -lgfortranbegin
#LIBPATHS = -L/usr/lib64 -L/afs/cern.ch/user/n/ndesai/public/lhapdf/lib/ \
#           -L/usr/libexec/CERNLIB/2006/lib

OBJECTS =  ddk23.o dmcarlo.o cuts.o kin.o ddk4.o Ctq6Pdf.o
PDFLIB = $(LIBLHA)

pphlnu.x: pphlnu.o ddk23.o dmcarlo.o cuts.o kin.o ddk4.o  Ctq6Pdf.o
	gfortran -o pphlnu.x $(OBJECTS) pphlnu.o $(LIBPATHS) $(DEBUG) 

pphlnu.o: pphlnu.f90
	gfortran -c pphlnu.f90 $(DEBUG)

higgsdecay.x:higgsdecay.o ddk23.o dmcarlo.o cuts.o kin.o ddk4.o
	gfortran -o higgsdecay.x  $(OBJECTS) higgsdecay.o

higgsdecay.o:higgsdecay.f90
	gfortran -c higgsdecay.f90 $(DEBUG)

zh.x: zh.o  dmcarlo.o ddk23.o cuts.o kin.o
	gfortran -o zh.x *.o 

zh.o: zh.f90
	gfortran -c zh.f90

muon.x: mudecay.o dmcarlo.o ddk23.o cuts.o kin.o
	gfortran -o muon.x *.o 

mudecay.o: mudecay.f90
	gfortran -c mudecay.f90 $(DEBUG)

Ctq6Pdf.o:  Ctq6Pdf.f
	gfortran -c  Ctq6Pdf.f $(DEBUG)

dmcarlo.o: dmcarlo.f
	gfortran -c dmcarlo.f $(DEBUG)

ddk4.o: ddk4.f
	gfortran -c ddk4.f $(DEBUG)

ddk23.o:  ddk23.f
	gfortran -c ddk23.f $(DEBUG)

cuts.o: cuts.f
	gfortran -c cuts.f $(DEBUG)

kin.o: kin.f
	gfortran -c kin.f $(DEBUG)

.PHONY:
clean:
	rm -f *.o
