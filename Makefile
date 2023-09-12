COMPILER=g++
LDFLAGS=-lcnpy -lz
STD=c++11
NAME=particles

all:
	$(COMPILER) $(NAME).cpp -o $(NAME) $(LDFLAGS) --std=$(STD)
