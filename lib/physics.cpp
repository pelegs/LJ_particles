#include <algorithm>
#include <cassert>
#include <cmath>
#include <glm/ext/matrix_transform.hpp>
#include <iostream>

// Physics-ish constants
const double R_CUTOFF = 5.0; // Cutoff distance for neighbor finding
const double GRAV = 1.0E2; // Gravitational constant
const double LJ_E = 1.0E6; // Lennard-Jones energy

// Functions
double distance1D(const double &x, const double &y) {
  // This is used for neighbor finding.
  // 1D distance is defined so that there's no need to use expensive squares and
  // square root in the case of distances in a single axis.
  return std::abs(x - y);
}

double U_LJ(double E, double S, double x) {
  // Lennard-Jones potential
  double S_x_6 = std::pow(S / x, 6.0);
  double S_x_12 = std::pow(S_x_6, 2.0);
  return 4.0 * E * (S_x_12 - S_x_6);
}

double F_LJ(double S, double x) {
  // Lennard-Jones force
  double S_6 = std::pow(S, 6.0);
  return -1.0 * 24.0 * LJ_E * S_6 * (std::pow(x, 6.0) - 2 * S_6) /
         std::pow(x, 13.0);
}

double F_HOOK(double K, double x, double x0) {
  // Hook's law for a simple, non-damped spring.
  return -1.0 * K * (x - x0);
}
