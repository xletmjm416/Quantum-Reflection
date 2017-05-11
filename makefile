CC=g++
CFLAGS=-I. -g
DEPS = matrix.h
OBJ = matrix.o 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

matrix: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)