COMPILER=g++
LDFLAGS=-lcnpy -lz
STD=c++11
MAIN=particles
LIBFOLDER=lib

physics:
	rm $(LIBFOLDER)/physics.o
	$(COMPILER) -c $(LIBFOLDER)/physics.cpp -o $(LIBFOLDER)/physics.o --std=$(STD)
	ar rvs $(LIBFOLDER)/physics.a $(LIBFOLDER)/physics.o

main:
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LIBFOLDER)/physics.a $(LDFLAGS) --std=$(STD)

debug:
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LDFLAGS) --std=$(STD) -g
