#include "spring.hpp"
#include "maths.hpp"
#include "physics.hpp"
#include "particles.hpp"
#include <iostream>

Spring::Spring(Particle &p1, Particle &p2, double K, double L) {
  this->p1 = &p1;
  this->p2 = &p2;
  this->K = K;
  if (L == -1.0)
    this->L = glm::distance(p1.get_pos(), p2.get_pos());
  else
    this->L = L;
}

void Spring::get_data() {
  std::cerr << this->K << " " << this->L << " " << this->p1->get_id() << " "
            << this->p2->get_id() << std::endl;
}

Particle *Spring::get_particle(int id) {
  if (id == 0)
    return this->p1;
  return p2;
}

double Spring::particles_distance() {
  return glm::distance(this->p1->get_pos(), this->p2->get_pos());
}

double Spring::hook_force(double x) {
  return F_HOOK(this->K, x, this->L);
}

void Spring::apply_force() {
  vec2 dr = this->p2->look_at(*this->p1);
  double x = this->particles_distance();
  vec2 F12 = F_HOOK(this->K, x, this->L) * dr;
  vec2 F21 = -1.0 * F12;
  this->p1->add_force(F12, MAX_FORCE);
  this->p2->add_force(F21, MAX_FORCE);
}
