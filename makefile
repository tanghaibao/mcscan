VER=0.8
PROG=mcscan
DIST=$(PROG)$(VER)
SRCS=struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc pog.cc permutation.cc
OBJS=$(SRCS:.cc=.o) 
CC=g++
CFLAGS=-O3 -Wall -ansi -pedantic-errors -I. -DVER=$(VER)

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

test:
	$(MAKE)
	python mcscan.py at_vv

dist:
	$(MAKE) clean
	cd .. && tar czf $(DIST).tar.gz $(DIST)/
