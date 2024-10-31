CC = g++
CFLAGS  = -O3 -Wall -Winline -Wshadow -fopenmp -std=c++17
TARGET = mgsolve

all: main.cpp 
	$(CC) $(CFLAGS) main.cpp -o $(TARGET)
