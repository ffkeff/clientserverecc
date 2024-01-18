CC=gcc
CFLAGS=-g -Wall -lgmp

all: client

ecc.o: ecc.h ecc.c
	$(CC) $(CFLAGS) -c ecc.c

client: ecc.o client.c
	$(CC) $(CFLAGS) -o client client.c ecc.o


