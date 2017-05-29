CC=g++
CFLAGS=-Wall -g -static-libstdc++
COMMON_FILES=common.h

main: main.o QMSystem.o FFT.o
	$(CC) $(CFLAGS) -o main main.o QMSystem.o FFT.o

main.o: main.cpp main.h $(COMMON_FILES)
	$(CC) $(CFLAGS) -c -o main.o main.cpp
	
QMSystem.o: QMSystem.cpp QMSystem.h $(COMMON_FILES)
	$(CC) $(CFLAGS) -c -o QMSystem.o QMSystem.cpp

FFT.o: FFT.cpp FFT.h $(COMMON_FILES)
	$(CC) $(CFLAGS) -c -o FFT.o FFT.cpp