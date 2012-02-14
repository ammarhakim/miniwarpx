include $(MINIWARPX)/src/make.in

INCLUDES += -I../core
TARGET = xminiwarpx

OBJS =  $(MINIWARPX)/src/lib/wave_source_advance.o \
	$(MINIWARPX)/src/lib/rkdg_limiter_characteristics.o \
	../core/rp3eu.o \
	../core/setprob.o \
	setkappa_1.o \
	transform_1.o \
	src.o \
	flux.o \
	bc.o \
	qinit.o

# standard library objects and target definitions
include $(MINIWARPX)/src/stdobjs.in

all : $(TARGET)
