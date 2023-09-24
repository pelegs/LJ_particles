#include "particles.hpp"
#include <vector>

class ParticleSystem {
  int num_particles;
  std::vector<Particle *> particle_list;

public:
  ParticleSystem();
  ~ParticleSystem();

  // Particle management
  void add_particle(Particle *p);
  void remove_particle();
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
  void append_new_data(std::vector<double> &data);
};
