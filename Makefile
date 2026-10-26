# makefile for alntools

CFLAGS = -O3 -DUTILS_THREADSAFE
#CFLAGS = -g -DUTILS_THREADSAFE	# for debugging

LIBS = -lpthread -lz -lm

ALL = tanbed gdbmask svfind taco ONEview tancons satmatch satscan

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

tancons.o: alntools.h ONElib.h $(UTILS_HEADERS)

gdb.o: alntools.h ONElib.h seqio.h $(UTILS_HEADERS)

satscan.o: onerep.h alntools.h satscore.h seqio.h ONElib.h $(UTILS_HEADERS)
satscore.o: satscore.h $(UTILS_HEADERS)

SEQIO_OPTS = -DONEIO
seqio.o: seqio.c seqio.h ONElib.h $(UTILS_HEADERS)
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

alnseq.o: alnseq.h ONElib.h

alncode.o: alncode.h align.h

### programs

tanbed: tanbed.o gdb.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

tancons: tancons.o gdb.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

gdbmask: gdb.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) -D GDB_MASK $(CFLAGS) -o $@ $^ $(LIBS)

taco: taco.c gdb.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

svfind: svfind.c alnseq.o alncode.o gdb.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

satmatch: satmatch.c satscore.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

satscan: satscan.o satscore.o gdb.o seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

### test

### end of file
