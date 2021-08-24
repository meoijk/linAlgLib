# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -Wall -g
 
# ****************************************************
# Targets needed to bring the executable up to date
 
TARGET = linAlgLib 

main: main.o 
	$(CC) $(CFLAGS) -o $(TARGET) main.o 
 
# The main.o target can be written more simply
 
main.o: main.cxx 
	$(CC) $(CFLAGS) -c main.cxx

clean:
	$(RM) $(TARGET) *.o 