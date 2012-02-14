include $(MINIWARPX)/src/make.in

INCLUDES += -I../core
TARGET = xminiwarpx

OBJS =  $(MINIWARPX)/src/lib/bc.o \
	$(MINIWARPX)/src/lib/wave_source_advance.o \
	$(MINIWARPX)/src/lib/setkappa.o \
	$(MINIWARPX)/src/lib/transform.o \
	../core/rkdg_twofluid_limiter.o \
	../core/rkdg_euler_limiter.o \
	../core/rkdg_maxwell_limiter.o \
	../core/rp_euler.o \
	../core/rp_maxwell.o \
	../core/rp_twofluid.o \
	../core/flux.o \
	../core/setprob.o \
	../core/src.o \
	qinit.o

#	../core/rkdg_twofluid_limiter.o \
#	../core/rkdg_euler_limiter.o \
#	../core/rkdg_maxwell_limiter.o \



# standard library objects and target definitions
include $(MINIWARPX)/src/stdobjs.in

all : $(TARGET)


