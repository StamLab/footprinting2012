CC	= g++
FLAGS	= -ansi -pedantic -Wall -static -O3 -s
DFLAGS	= -ansi -pedantic -Wall -static -O0 -g

NAME1	= fp2012
SOURCE1	= $(NAME1).cpp
BINDIR	= ../bin

detect:
	mkdir -p $(BINDIR) && $(CC) $(FLAGS) $(SOURCE1) -o $(BINDIR)/$(NAME1)

debug:
	mkdir -p $(BINDIR) && $(CC) $(DFLAGS) $(SOURCE1) -o $(BINDIR)/debug.$(NAME1)

clean:
	rm -f $(BINDIR)/$(NAME1)
	rm -f $(BINDIR)/debug.$(NAME1)
