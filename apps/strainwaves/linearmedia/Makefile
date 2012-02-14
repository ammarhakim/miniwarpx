include $(MINIWARPX)/src/make.in

INCLUDES += -I../core
TARGET = xminiwarpx

OBJS =  $(MINIWARPX)/src/lib/src.o \
	$(MINIWARPX)/src/lib/wave_source_advance.o \
	$(MINIWARPX)/src/lib/setkappa.o \
	$(MINIWARPX)/src/lib/transform.o \
	$(MINIWARPX)/src/lib/after_step.o \
	$(MINIWARPX)/src/lib/rkdg_limiter_characteristics.o \
	../core/rp.o \
	../core/flux.o \
	../core/setprob.o \
	bc.o \
	qinit.o

# standard library objects and target definitions
include $(MINIWARPX)/src/stdobjs.in

all : $(TARGET)


