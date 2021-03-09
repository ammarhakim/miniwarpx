# Amalgamated Makefile, generated automatically. Set compiler or
# modify as needed. For example, to set compiler do:
#
# make CC=mpicc 
#

MINIWARPX = ${HOME}/research/ancient-history/miniwarpx
CFLAGS = -O3 -g -I. -I${MINIWARPX}/src/etc -I${MINIWARPX}/src/lib
CXXFLAGS = -O3 -g -I. -I${MINIWARPX}/src/etc -I${MINIWARPX}/src/lib
LINK_DIR = -L${MINIWARPX}/src

all: xwarpx

xwarpxobjs = qinit.o ../core/flux.o ../core/rp.o

xwarpx: ${xwarpxobjs}
	${CXX} -O3 -g qinit.o ../core/flux.o ../core/rp.o ${LINK_DIR} -lwarpx -o xwarpx

clean:
	rm -rf xwarpx ${xwarpxobjs}
