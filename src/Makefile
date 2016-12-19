OBJS  =	main.o
SOURCE  =	main.cpp
HEADER  =	IO.h	memory.h	hash.h  hypercube.h
OUT   =	dolphinn
CXX =	g++
FLAGS	=	-pthread    -std=c++0x	-Wall   -O3 -Qunused-arguments

all:	$(OBJS)
	$(CXX)	$(OBJS)	-o	$(OUT)	$(FLAGS)
	make	-f	Makefile	clean
  
# create/compile the individual files >>separately<< 
main.o:	main.cpp
	$(CXX)	-c	main.cpp	$(FLAGS)
    
.PHONY:	all
# clean house
clean:
	rm -f $(OBJS)

# do a bit of accounting
count:
	wc $(SOURCE) $(HEADER)
