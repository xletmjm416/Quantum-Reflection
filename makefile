CC=g++
CFLAGS=-Wall -g -static-libstdc++

main: main.o QMSystem.o
	$(CC) $(CFLAGS) -o main main.o QMSystem.o
	
QMSystem.o: QMSystem.cpp QMSystem.h
	$(CC) $(CFLAGS) -c -o QMSystem.o QMSystem.cpp
	
main.o: main.cpp QMSystem.h
	$(CC) $(CFLAGS) -c -o main.o main.cpp