CXX = gcc
EXEC = exe
CXXFLAGS += -pedantic -Wall -MMD -MP
LFLAGS += -lgsl -lgslcblas -lm

SRC=$(wildcard *.c)
OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)

$(EXEC): $(OBJ) 
	$(CXX) $^ -o $@ $(LFLAGS)

clean:
	rm -f $(OBJ) $(DEP)

run: 
	./$(EXEC)
	
.PHONY: clean run
