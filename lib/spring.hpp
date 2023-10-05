#ifndef SPRING
#define SPRING

#include <iostream>
#include "maths.hpp"
#include "physics.hpp"
#include "particles.hpp"
 
class Spring {
  Particle *p1, *p2;
  double K, L;

public:
  Spring(Particle &p1, Particle &p2, double K, double L);
  void get_data();
  Particle *get_particle(int id);
  double particles_distance();
  double hook_force(double x);
  void apply_force();
};

#endif // !SPRING
