include $(MINIWARPX)/src/make.in

INCLUDES += -I../core
TARGET = xminiwarpx

OBJS =  $(MINIWARPX)/src/lib/bc.o \
	$(MINIWARPX)/src/lib/wave_source_advance.o \
	$(MINIWARPX)/src/lib/src.o \
	$(MINIWARPX)/src/lib/rkdg_limiter_characteristics.o \
	setkappa.o \
	transform.o \
	../core/rp3eu.o \
	../core/flux.o \
	../core/setprob.o \
	qinit.o

# standard library objects and target definitions
include $(MINIWARPX)/src/stdobjs.in

all : $(TARGET)


