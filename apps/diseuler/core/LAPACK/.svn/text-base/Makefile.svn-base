CXX = g++
CC = gcc
LD = ld

## uncomment appropriate FLAG to optimize/debug
## optmize:
FLAGS = -Wall -O2
#FLAGS = -g -Wall

LIB_LAPACK_OBJS = blas_lite.o dlapack_lite.o f2c_lite.o zlapack_lite.o

LINK_DIRS = -L.
LINK_FLAGS = -llapack_lite

.c.o:
	$(CXX) -c $(INCLUDES) $(FLAGS) -o $@ $<

.cc.o:
	$(CXX) -c $(INCLUDES) $(FLAGS) -o $@ $<
	$(CXX)  $(FLAGS) $(OBJS) $(LINK_DIRS) $(LINK_FLAGS) -o $(TARGET)

lapack: $(LIB_LAPACK_OBJS)

clean:
	-rm -f *.o *.so *.a