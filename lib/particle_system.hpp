#include "cnpy.h"
#include "particles.hpp"
#include <vector>

class ParticleSystem {
  unsigned long num_particles, num_steps;
  vec2 space_dimensions;
  std::vector<Particle *> particle_list = {};
  std::vector<Particle *> particle_list_sorted[2];
  std::vector<Wall *> walls = {};
  std::vector<double> trajectories = {}, AABB_min = {}, AABB_max = {}, forces = {};
  std::vector<int> neighbors_matrix = {};
  std::vector<int> sorted_particle_ids_X, sorted_particle_ids_Y; // temp

public:
  ParticleSystem();
  ParticleSystem(double width, double height);
  ~ParticleSystem();

  // Particle management
  void add_particle(Particle *p);
  // void remove_particle();
  Particle *get_particle(int i);
  std::vector<Particle *> get_particle_list();
  std::vector<Wall *> get_wall_list();

  // Wall managment
  void add_wall(Wall *wall);

  // Collision detection
  void reset_neighbors();
  void sort_particles_by_min_AABB(int axis);
  void sort_particles_all_directions();
  void assign_neighbors_by_axis(int axis);
  void assign_neighbors();
  void interact_with_walls(double atol);

  // Dynamics
  void calc_new_positions(const double &dt, const bool &update_data, const bool &update_trajectories_data);
  void calc_accelerations();
  void calc_new_velocities(const double &dt);
  void move_particles(const double &dt, const bool &with_interactions);
  void interact(bool LJ, bool gravity, bool springs);

  // Data managment
  void update_neighbors_matrix();
  void validate_neighbors();
  void save_data(std::string filename, bool save_particle_data,
                 bool save_neighbor_matrix, bool save_sort_data, bool save_AABB_data, bool save_forces);
};
