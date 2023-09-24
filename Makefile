COMPILER=g++
LDFLAGS=-lcnpy -lz
STD=c++11
MAIN=main
LIBFOLDER=lib

maths:
	$(COMPILER) -c $(LIBFOLDER)/maths.cpp -o $(LIBFOLDER)/maths.o --std=$(STD)

physics: maths
	$(COMPILER) -c $(LIBFOLDER)/physics.cpp -o $(LIBFOLDER)/physics.o --std=$(STD)

otherfuncs:
	$(COMPILER) -c $(LIBFOLDER)/otherfuncs.cpp -o $(LIBFOLDER)/otherfuncs.o --std=$(STD)

particles: physics otherfuncs
	$(COMPILER) -c $(LIBFOLDER)/particles.cpp -o $(LIBFOLDER)/particles.o --std=$(STD)

main: physics otherfuncs particles
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LIBFOLDER)/maths.o $(LIBFOLDER)/physics.o $(LDFLAGS) --std=$(STD)

debug:
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LDFLAGS) --std=$(STD) -g
