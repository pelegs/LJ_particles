#ifndef PHYSICS
#define PHYSICS

// Physics-ish constants
const double R_CUTOFF = 10.0; // Cutoff distance for neighbor finding
const double GRAV = 1.0E2; // Gravitational constant
const double LJ_E = 1.0E6; // Lennard-Jones energy

// Funcs
double distance1D(const double &x, const double &y);
double U_LJ(double E, double S, double x);
double F_LJ(double S, double x);
double F_HOOK(double K, double x, double x0);

#endif // !PHYSICS
