include $(MINIWARPX)/src/make.in

OBJS =	flux.o \
	setprob.o \
	rp_euler.o \
	rp_maxwell.o \
	rp_twofluid.o \
	rkdg_twofluid_limiter.o \
	rkdg_euler_limiter.o \
	rkdg_maxwell_limiter.o \
	src.o

all : $(OBJS)

clean:
	-rm -f *.o
count:
	wc -l *.cc *.h
