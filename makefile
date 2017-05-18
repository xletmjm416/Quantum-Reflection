CC=g++
CFLAGS=-g -static-libstdc++
DEPS = 
OBJ = base.o 

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

base: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)