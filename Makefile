#compile

CC = g++
LD = g++
AR = ar


#linkedit
DEBUG = -D _MY_DEBUG_
CFLAGS = -Wall -g
LDFLAGS = -Wall

AFLAGS = -r

SOURCES = main.cpp utilities.cpp basicMatrix.cpp linearEquation.cpp mmse.cpp nonLinearSolutions.cpp
LISTA1MAIN = lista1Main
LISTAMAINOBJS = $(SOURCES:.c=.o)

EXECS = $(LISTA1MAIN)

ALL = $(LISTA1MAIN)

$(LISTA1MAIN) : $(LISTAMAINOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o $(ALL)

