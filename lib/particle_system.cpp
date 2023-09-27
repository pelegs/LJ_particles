#include "particle_system.hpp"
#include <algorithm>
#include <iostream>

ParticleSystem::ParticleSystem() {
  this->num_particles = 0;
  this->num_steps = 0;
  this->particle_list.clear();
}

ParticleSystem::~ParticleSystem() {
  for (auto particle : this->particle_list)
    delete particle;
}

// Particle management
void ParticleSystem::add_particle(Particle *p) {
  this->particle_list.push_back(p);
  this->particle_list_sorted[X].push_back(p);
  this->particle_list_sorted[Y].push_back(p);
  this->num_particles++;
}

// void ParticleSystem::remove_particle() {
//   this->num_particles--;
// } // TBW

Particle *ParticleSystem::get_particle(int i) { return this->particle_list[i]; }

std::vector<Particle *> ParticleSystem::get_particle_list() {
  return this->particle_list;
}

// Collision detection
void ParticleSystem::reset_neighbors() {
  for (auto particle : this->particle_list)
    particle->reset_neighbors();
}

void ParticleSystem::sort_particles(int axis) {
  std::sort(this->particle_list_sorted[axis].begin(),
            this->particle_list_sorted[axis].end(), CompareParticlesAABB(axis));
}

void ParticleSystem::sort_particles_all_directions() {
  this->sort_particles(X);
  this->sort_particles(Y);

  // temp
  this->sorted_particle_ids_X.clear();
  this->sorted_particle_ids_Y.clear();
  for (auto particle : this->particle_list_sorted[X])
    this->sorted_particle_ids_X.push_back(particle->get_id());
  for (auto particle : this->particle_list_sorted[Y])
    this->sorted_particle_ids_Y.push_back(particle->get_id());
}

void ParticleSystem::assign_neighbors_by_axis(int axis) {
  for (auto it1 = this->particle_list_sorted[axis].begin();
       it1 != std::prev(this->particle_list_sorted[axis].end()); ++it1) {
    auto &p1 = *it1;
    for (auto it2 = (it1 + 1);
         it2 != std::prev(this->particle_list_sorted[axis].end()); ++it2) {
      auto &p2 = *it2;
      if (p2->get_min_AABB(axis) < p1->get_max_AABB(axis)) {
        p1->add_neighbor(axis, p2);
        p2->add_neighbor(axis, p1);
      } else {
        break;
      }
    }
  }
}

void ParticleSystem::assign_neighbors() {
  this->assign_neighbors_by_axis(X);
  this->assign_neighbors_by_axis(Y);
  for (auto particle : this->particle_list) {
    particle->generate_neighbors_list_by_intersection();
  }
}

// Dynamics
void ParticleSystem::calc_new_positions(const double &dt) {
  for (auto &p : this->particle_list) {
    p->calc_new_pos(dt);
    this->trajectories.push_back(p->get_x());
    this->trajectories.push_back(p->get_y());
  }
  this->num_steps++;
}

void ParticleSystem::calc_accelerations() {
  for (auto &p : this->particle_list) {
    p->calc_acc();
  }
}

void ParticleSystem::calc_new_velocities(const double &dt) {
  for (auto &p : this->particle_list) {
    p->calc_new_vel(dt);
    p->check_wall_collision(this->space_dimensions[X],
                            this->space_dimensions[Y]);
  }
}

void ParticleSystem::move_particles(const double &dt) {
  calc_new_positions(dt);
  // interaction
  calc_accelerations();
  calc_new_velocities(dt);
}

// Data managment
void ParticleSystem::update_trajectory_data() {
  for (auto particle : this->particle_list) {
    this->trajectories.push_back(particle->get_x());
    this->trajectories.push_back(particle->get_y());
  }
}

void ParticleSystem::update_neighbors_matrix() {
  std::vector<int> neighbors_row(this->num_particles, 0);
  for (auto particle : this->particle_list) {
    std::fill(neighbors_row.begin(), neighbors_row.end(), 0);
    for (auto neighbor : particle->get_neighbors_list())
      neighbors_row[neighbor->get_id()] = 1;
    this->neighbors_matrix.insert(this->neighbors_matrix.begin(),
                                  neighbors_row.begin(), neighbors_row.end());
  }
}

void ParticleSystem::save_data(std::string filename,
                               bool save_particle_data = false,
                               bool save_neighbor_matrix = false,
                               bool save_sort_data = false) {
  cnpy::npz_save(filename, "space_dimensions", &this->space_dimensions[0], {2},
                 "w");

  std::vector<unsigned long> params = {this->num_particles, this->num_steps};
  cnpy::npz_save(filename, "parameters", &params[0], {2}, "a");

  if (this->trajectories.size()) {
    cnpy::npz_save(filename, "trajectories", &this->trajectories[0],
                   {this->num_steps, this->num_particles, 2}, "a");
  }

  if (this->trajectories.size()) {
    cnpy::npz_save(filename, "trajectories", &this->trajectories[0],
                   {this->num_steps, this->num_particles, 2}, "a");
  }

  if (save_particle_data) {
    std::vector<double> masses;
    std::vector<double> radii;
    std::vector<double> bounding_distances;
    for (auto particle : this->particle_list) {
      masses.push_back(particle->get_mass());
      radii.push_back(particle->get_radius());
      bounding_distances.push_back(particle->get_bounding_distance());
    }
    cnpy::npz_save(filename, "masses", &masses[0], {this->num_particles}, "a");
    cnpy::npz_save(filename, "radii", &radii[0], {this->num_particles}, "a");
    cnpy::npz_save(filename, "bounding_distances", &bounding_distances[0],
                   {this->num_particles}, "a");
  }

  if (save_neighbor_matrix && this->neighbors_matrix.size()) {
    cnpy::npz_save(filename, "neighbors_matrix", &this->neighbors_matrix[0],
                   {this->num_steps, this->num_particles, this->num_particles},
                   "a");
  }

  if (save_sort_data) {
    cnpy::npz_save(filename, "sort_by_x", &this->sorted_particle_ids_X[0],
                   {this->num_particles}, "a");
    cnpy::npz_save(filename, "sort_by_y", &this->sorted_particle_ids_Y[0],
                   {this->num_particles}, "a");
  }
}
