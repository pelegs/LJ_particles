COMPILER=g++
LDFLAGS=-lcnpy -lz
STD=c++11
MAIN=main
LIBFOLDER=lib

maths:
	$(COMPILER) -c $(LIBFOLDER)/maths.cpp -o $(LIBFOLDER)/maths.o --std=$(STD)

physics: maths
	$(COMPILER) -c $(LIBFOLDER)/physics.cpp -o $(LIBFOLDER)/physics.o --std=$(STD)

main: physics
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LIBFOLDER)/maths.o $(LIBFOLDER)/physics.o $(LDFLAGS) --std=$(STD)

debug:
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LDFLAGS) --std=$(STD) -g
