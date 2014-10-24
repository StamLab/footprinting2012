CC	= g++
FLAGS	= -ansi -pedantic -Wall -static -O3 -s
DFLAGS	= -ansi -pedantic -Wall -static -O0 -g

SOURCE1	= detect-cache.cpp
BINDIR	= ../bin

detect:
	mkdir -p $(BINDIR) && $(CC) $(FLAGS) $(SOURCE1) -o $(BINDIR)/detect-cache

debug:
	mkdir -p $(BINDIR) && $(CC) $(DFLAGS) $(SOURCE1) -o $(BINDIR)/detect-cache-debug

clean:
	rm -f $(BINDIR)/detect-cache
