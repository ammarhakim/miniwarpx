include $(MINIWARPX)/src/make.in

HEADERS = section.h
OBJS = 	getline.o \
	section.o \
	inpparse.o \
	extractors.o \
	copy_file.o

all : $(OBJS)

$(TARGET): $(OBJS) $(HEADERS)
	$(CC)  $(FLAGS) $(OBJS) $(LINK_DIRS) $(LINK_FLAGS) -o $(TARGET)

clean:
	-rm -f *.o

count:
	-wc -l *.c *.cc *.h
