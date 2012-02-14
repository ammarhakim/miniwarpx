include $(MINIWARPX)/src/make.in

INCLUDES += -I../core -I../core/LAPACK
TARGET = xminiwarpx

OBJS =  $(MINIWARPX)/src/lib/bc.o \
	$(MINIWARPX)/src/lib/wave_source_advance.o \
	$(MINIWARPX)/src/lib/setkappa.o \
	$(MINIWARPX)/src/lib/transform.o \
	$(MINIWARPX)/src/lib/after_step.o \
	../core/rkdg_euler_limiter.o \
	../core/rp3eu.o \
	../core/flux.o \
	../core/setprob.o \
	../core/src.o \
	../core/LAPACK/dlapack_lite.o \
	../core/LAPACK/blas_lite.o \
	../core/LAPACK/f2c_lite.o \
	qinit.o

# standard library objects and target definitions
include $(MINIWARPX)/src/stdobjs.in

all : $(TARGET)


