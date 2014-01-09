## CS130 Fall 2013: Computer Graphics
## Lab 1 Makefile
##
## This file does not need to be modified
#########################################

LIBS = -lglut -lGLU -lGL -lXmu -lXext -lXi -lX11 -lm
CC = g++

## Global header files
INCLUDE =

## Object files and executables
MAIN_OUT = render

## Requirements for each command
MAIN_REQS = proj2.cpp

## Targets to compile for each command
MAIN_TARGETS = proj2.cpp

all: main

## Main 
main: $(MAIN_REQS)
	$(CC) $(MAIN_TARGETS) $(LIBS) -o $(MAIN_OUT)

clean:
	rm -f *~ *.o *.out
