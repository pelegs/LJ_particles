#include "particle_system.hpp"

ParticleSystem::ParticleSystem() {
  this->num_particles = 0;
  this->particle_list.clear();
}

ParticleSystem::~ParticleSystem() {
  for (auto particle : this->particle_list)
    delete particle;
}

// Particle management
void ParticleSystem::add_particle(Particle *p) {
  this->particle_list.push_back(p);
  this->num_particles++;
}

void ParticleSystem::remove_particle() {
  this->num_particles--;
} // TBW

Particle* ParticleSystem::get_particle(int i) { return this->particle_list[i]; }

std::vector<Particle *> ParticleSystem::get_particle_list() {
  return this->particle_list;
}

// Dynamics
void ParticleSystem::calc_new_positions(const double &dt) {
  for (auto &p : this->particle_list) {
    p->calc_new_pos(dt);
  }
}

void ParticleSystem::calc_accelerations() {
  for (auto &p : this->particle_list) {
    p->calc_acc();
  }
}

void ParticleSystem::calc_new_velocities(const double &dt, const double &width,
                                         const double &height) {
  for (auto &p : this->particle_list) {
    p->calc_new_vel(dt);
    p->check_wall_collision(width, height);
  }
}

void ParticleSystem::move_particles(const double &dt, const double &width,
                                    const double &height) {
  calc_new_positions(dt);
  // interaction
  calc_accelerations();
  calc_new_velocities(dt, width, height);
}

// Data managment
void ParticleSystem::save_trajectory_data() {
  for (auto particle : this->particle_list) {
    this->trajectories.push_back(particle->get_x());
    this->trajectories.push_back(particle->get_y());
  }
}
