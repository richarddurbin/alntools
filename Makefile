# makefile for alntools

CFLAGS = -O3
#CFLAGS = -g	# for debugging

LIBS = -lpthread -lz

ALL = tanbed ONEview

DESTDIR = ~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	$(RM) -r *.dSYM

### object files

UTILS_OBJS = utils.o array.o dict.o
UTILS_HEADERS = utils.h array.h dict.h
$(UTILS_OBJS): $(UTILS_HEADERS)

ONElib.o: ONElib.h 

tanbed.o: alntools.h ONElib.h $(UTILS_HEADERS)

### programs

tanbed: tanbed.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

### test

### end of file
