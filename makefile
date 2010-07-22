VER=0.8
PROG=mcscan
DIST=$(PROG)-$(VER)
SRCS=basic.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc pog.cc permutation.cc
OBJS=$(SRCS:.cc=.o) 
CC=g++
CFLAGS=-O3 -Wall -ansi -pedantic-errors -I. -DVER=$(VER)
CFLAGS+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(PROG)
	rm -f *.o
	rm -f data/*.aligns
	rm -f data/*.blocks

test: $(PROG)
	run.sh

dist:
	$(MAKE) clean
	git archive HEAD | gzip > ../$(DIST).tar.gz

