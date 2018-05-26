# Makefile for sierpinski gasket program

CC = gcc
CFLAGS = -ansi -pedantic -Werror -Wall -O2 -D_GNU_SOURCE
LFLAGS = -lm -ltiff -lpthread

all: sierpinski

sierpinski: sierpinski_gasket.c sierpinski.h
	$(CC) $(CFLAGS) sierpinski_gasket.c -o sierpinski $(LFLAGS)

clean:
	\rm sierpinski
