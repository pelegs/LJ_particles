COMPILER=g++
LDFLAGS=-lcnpy -lz
STD=c++17
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

particle_system: particles
	$(COMPILER) -c $(LIBFOLDER)/particle_system.cpp -o $(LIBFOLDER)/particle_system.o --std=$(STD)

spring: particles
	$(COMPILER) -c $(LIBFOLDER)/spring.cpp -o $(LIBFOLDER)/spring.o --std=$(STD)

main: physics otherfuncs particle_system spring
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LIBFOLDER)/maths.o $(LIBFOLDER)/physics.o $(LIBFOLDER)/otherfuncs.o $(LIBFOLDER)/particles.o $(LIBFOLDER)/particle_system.o $(LIBFOLDER)/spring.o $(LDFLAGS) --std=$(STD)

debug:
	$(COMPILER) $(MAIN).cpp -o $(MAIN) $(LDFLAGS) --std=$(STD) -g
