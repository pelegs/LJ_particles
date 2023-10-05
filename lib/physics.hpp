#ifndef PHYSICS
#define PHYSICS

#include "maths.hpp"
#include <glm/ext/matrix_transform.hpp>

const double LJ_E = 1.0E6;
const double GRAV = 1.0E5;
const double MAX_FORCE = 1.0E3;

double U_LJ(double, double, double);
double F_LJ(double, double);
double F_HOOK(double, double, double);

#endif // !PHYSICS
