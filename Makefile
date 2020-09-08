#compile

CC = g++
LD = g++
AR = ar


#linkedit
DEBUG = -D _MY_DEBUG_
CFLAGS = -Wall -O3
LFLAGS = -Wall

AFLAGS = -r


LISTAMAIN = lista1Main
LISTAMAINOBJS = basicMatrix.o main.o linearEquation.o utilities.o

EXECS = $(LISTAMAIN)

ALL = $(LISTAMAIN)

#Regra Implicita
.c.o: $(CC) $(CFLAGS) -c $<

lista1Main : $(LISTAMAINOBJS)
	$(LD) $(LFLAGS) -o $@ $(LISTAMAINOBJS)

clean:
	rm -f *.o $(ALL)

