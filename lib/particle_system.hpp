#include "particles.hpp"
#include <vector>

class ParticleSystem {
  int num_particles, num_steps;
  vec2 size;
  std::vector<Particle *> particle_list = {};
  std::vector<double> trajectories = {};

public:
  ParticleSystem();
  ~ParticleSystem();

  // Particle management
  void add_particle(Particle *p);
  // void remove_particle();
  Particle* get_particle(int i);
  std::vector<Particle *> get_particle_list();

  // Dynamics
  void calc_new_positions(const double &dt);
  void calc_accelerations();
  void calc_new_velocities(const double &dt, const double &width,
                           const double &height);
  void move_particles(const double &dt, const double &width,
                      const double &height);

  // Data managment
  void update_trajectory_data();
  void save_data(std::string filename);
};
