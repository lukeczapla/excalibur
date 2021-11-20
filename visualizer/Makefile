# Visualization software Makefile

CC=g++
PLATFORM=linux 		# either "osx" or "linux"

VIS_EXECUTABLE=visualizer

CFLAGS=-O3 -Wall

LINUXLIB=-lglut -lGL -lGLU -lm
OSXLIB=-framework opengl -framework glut -lm


all: vis_$(PLATFORM)

osx: vis_osx
linux: vis_linux

clean: 
	rm $(VIS_EXECUTABLE)

test: 
	$(VIS_EXECUTABLE) 7001


# *** Source code compilation

VIS_SOURCES=vis/vis.cpp vis/gl.cpp util/misc.cpp util/thread.cpp util/socket.cpp

vis_osx: $(VIS_SOURCES)
	$(CC) $(CFLAGS) $(VIS_SOURCES) -o $(VIS_EXECUTABLE) $(OSXLIB)

vis_linux: $(VIS_SOURCES)
	$(CC) $(CFLAGS) $(VIS_SOURCES) -o $(VIS_EXECUTABLE) $(LINUXLIB)


