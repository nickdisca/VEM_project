TARGET = main_project
CC = g++
STANDARD = -std=c++11
CFLAGS = -O2 #-Wall
INCLUDE = -I. -I/usr/local/include/eigen3

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.hpp)
SOURCES = $(wildcard *.cpp)

%.o: %.cpp $(HEADERS)
	$(CC) -c $(STANDARD) $(INCLUDE) $(CFLAGS) $*.cpp

$(TARGET): $(OBJECTS)
	$(CC) $(STANDARD) $(INCLUDE) $(OBJECTS) -o $(TARGET)

clean:
	-rm -f *.o
	-rm -f $(TARGET)
