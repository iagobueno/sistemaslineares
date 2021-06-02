#GRR20190171 Carlos Iago Bueno
CC = gcc
CFLAGS = -I
DEPS = functions.h
OBJ = labSisLin.o SistemasLineares.o utils.o functions.o
LIBS = -lm

labSisLin: $(OBJ)
	$(CC) $(OBJ) -o labSisLin $(LIBS)

labSisLin.o: labSisLin.c
	$(CC) -c labSisLin.c

utils.o: utils.c
	$(CC) -c utils.c

functions.o: functions.c
	$(CC) -c functions.c

SistemasLineares.o: SistemasLineares.c functions.o
	$(CC) -c SistemasLineares.c

clean:
	-rm -f $(OBJ)

purge: clean
	-rm -f labSisLin

ex: labSisLin
	-./labSisLin < teste

j: labSisLin
	-./labSisLin < jacob

debug: labSisLin
	$(CC) -g $(OBJ) -o debug

gdb: debug
	gdb ./labSisLin