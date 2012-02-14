include $(MINIWARPX)/src/make.in

HEADERS = farray.h 
OBJS = 	farray.o \
	utils.o \
	gauleg.o \
	out.o \
	write_grid.o \
	solve.o \
	advance.o \
	eval_fluctuations_qwave.o \
	eval_fluctuations_fwave.o \
	eval_fluctuations.o \
	wave_initialize.o \
	wave_step.o \
	wave_setup.o \
	wave_source_advance.o \
	wave_limiter.o \
	rkdg_initialize.o \
	rkdg_step.o \
	rkdg_setup.o \
	rkdg_algo.o \
	rkdg_eval_expansion.o \
	rkdg_limiter.o \
	rkdg_limiter_characteristics.o \
	rp.o \
	qinit.o \
	flux.o \
	src.o \
	rp.o \
	bc.o \
	read_inp.o \
	driver.o \
	setprob.o \
	setkappa.o \
	transform.o \
	init_output.o \
	before_step.o \
	after_step.o \
	main.o

all : $(OBJS)

$(TARGET): $(OBJS) $(HEADERS)
	$(CXX)  $(FLAGS) $(OBJS) $(LINK_DIRS) $(LINK_FLAGS) -o $(TARGET)

clean:
	-rm -f *.o
count:
	-wc -l *.cc *.h
