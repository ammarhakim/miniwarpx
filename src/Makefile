# Amalgamated Makefile, generated automatically. Set compiler or
# modify as needed. For example, to set compiler do:
#
# make CC=mpicc 
#

CFLAGS = -O3 -g -I. -Ietc -Ilib
CXXFLAGS = -O3 -g -I. -Ietc -Ilib

headers = etc/extractors.h etc/inpparse.h etc/section.h lib/farray.h lib/maccor_algo.h lib/miniwarpx.h lib/rkdg_algo.h lib/utils.h lib/wave_algo.h lib/wxarray.h lib/wxindexer.h lib/wxrange.h lib/wxsequencer.h

libobjs = etc/copy_file.o etc/extractors.o etc/getline.o etc/inpparse.o etc/section.o lib/advance.o lib/after_step.o lib/bc.o lib/before_step.o lib/driver.o lib/eval_fluctuations.o lib/eval_fluctuations_fwave.o lib/eval_fluctuations_qwave.o lib/flux.o lib/gauleg.o lib/init_output.o lib/maccor2_initialize.o lib/maccor2_setup.o lib/maccor2_source_advance.o lib/maccor2_step.o lib/main.o lib/maxs.o lib/out.o lib/out_dot_gkyl.o lib/qinit.o lib/read_inp.o lib/rkdg_algo.o lib/rkdg_eval_expansion.o lib/rkdg_initialize.o lib/rkdg_limiter.o lib/rkdg_limiter_characteristics.o lib/rkdg_setup.o lib/rkdg_step.o lib/rp.o lib/setkappa.o lib/setprob.o lib/solve.o lib/src.o lib/transform.o lib/utils.o lib/wave_initialize.o lib/wave_limiter.o lib/wave_setup.o lib/wave_source_advance.o lib/wave_step.o lib/write_grid.o lib/wxindexer.o lib/wxrange.o lib/wxsequencer.o

all: libwarpx.a

libwarpx.a: ${libobjs}
	ar -crs libwarpx.a ${libobjs}

clean:
	rm -rf libwarpx.a ${libobjs}
