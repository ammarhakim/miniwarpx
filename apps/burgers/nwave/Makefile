include $(MINIWARPX)/src/make.in

INCLUDES += -I../core
TARGET = xminiwarpx

OBJS =  $(MINIWARPX)/src/lib/bc.o \
	$(MINIWARPX)/src/lib/wave_source_advance.o \
	$(MINIWARPX)/src/lib/src.o \
	$(MINIWARPX)/src/lib/setprob.o \
	$(MINIWARPX)/src/lib/rkdg_limiter_characteristics.o \
	../core/rp.o \
	../core/flux.o \
	qinit.o

# standard library objects and target definitions
include $(MINIWARPX)/src/stdobjs.in

all : $(TARGET)


