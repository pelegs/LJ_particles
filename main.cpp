#include "cnpy.h"
#include <cassert>
#include <glm/ext/matrix_transform.hpp>
#include <indicators/progress_bar.hpp>
#include <iterator>
#include <random>

// Own lib
#include "lib/maths.hpp"
#include "lib/physics.hpp"
#include "lib/otherfuncs.hpp"
#include "lib/particles.hpp"
#include "lib/particle_system.hpp"
#include "lib/spring.hpp"

// GLM-related
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#define assertm(exp, msg)                                                      \
  assert(((void)msg, exp)) // use (void) to silence unused warnings

/******************************************/
/*        Class-relevant functions        */
/******************************************/

bool comapreParticleByXPos(const Particle *lhs, const Particle *rhs) {
  return lhs->get_x() < rhs->get_x();
}

bool comapreParticleByYPos(const Particle *lhs, const Particle *rhs) {
  return lhs->get_y() < rhs->get_y();
}

/**********************/
/*        Main        */
/**********************/

int main(int argc, char *argv[]) {
  std::cout << R_CUTOFF << std::endl;
  // randomness!
  srand(time(NULL));

  // params
  double width = atof(argv[1]);
  double height = atof(argv[2]);
  int num_particles = atoi(argv[3]);
  int num_steps = atoi(argv[4]);
  int skip = atoi(argv[5]);
  double dt = atof(argv[6]);
  std::string filename = argv[7];
  filename += ".npz";

  // Neighbor data
  std::vector<int> neighbors_matrix = {};
  std::vector<int> nmat_row(num_particles, 0);

  // init particles
  ParticleSystem particle_system;
  for (int i = 0; i < num_particles; i++) {
    particle_system.add_particle(new Particle(i, O_, O_, 1.0, 3.0));
  }
  particle_system.get_particle(20)->set_vel(100.0, 150.0);
  // Set special big particle
  // particle_system.get_particle(49)->set_vel(.0, .0);
  // particle_system.get_particle(49)->set_mass(3.0);
  // particle_system.get_particle(49)->set_radius(5.0);

  // Create mass and radius vectors for saving data
  std::vector<double> masses = {};
  std::vector<double> radii = {};
  for (auto particle : particle_system.get_particle_list()) {
    masses.push_back(particle->get_mass());
    radii.push_back(particle->get_radius());
  }

  // Create springs?
  int bounded_indices[4] = {20, 21, 28, 29};
  std::vector<Spring *> springs = {};
  for (int id1 : bounded_indices)
    for (int id2 : bounded_indices)
      if (id1 != id2) {
        springs.push_back(new Spring(*particle_system.get_particle(id1),
                                     *particle_system.get_particle(id2), 30.0,
                                     -1.0));
      }

  // springs.push_back(new Spring(*particle_system.get_particle(id),
  //                              *particle_system.get_particle(21), 50.0,
  //                              -1.0));

  // Sorted vecs
  std::vector<Particle *> particles_x_pos = {};
  std::vector<Particle *> particles_y_pos = {};
  for (auto particle : particle_system.get_particle_list()) {
    particles_x_pos.push_back(particle);
    particles_y_pos.push_back(particle);
  }

  // Progress bar
  int progress_perc = 0;
  std::string pgtext = "Simulation running: 0/" + std::to_string(num_particles);
  indicators::ProgressBar bar{
      indicators::option::BarWidth{100},
      indicators::option::Start{"["},
      indicators::option::Fill{"■"},
      indicators::option::Lead{"■"},
      indicators::option::Remainder{"-"},
      indicators::option::End{" ]"},
      indicators::option::PostfixText{pgtext},
      indicators::option::FontStyles{
          std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};

  // Simulation
  std::vector<double> trajectories = {};
  std::vector<double> x_pos = {};
  std::vector<double> y_pos = {};

  // Neighor finding
  for (int step = 0; step < num_steps; step++) {
    // Sort x- and y-position vectors
    std::sort(particles_x_pos.begin(), particles_x_pos.end(),
              comapreParticleByXPos);
    std::sort(particles_y_pos.begin(), particles_y_pos.end(),
              comapreParticleByYPos);

    // Create neighbor lists for particles
    for (auto particle : particle_system.get_particle_list())
      particle->reset_neighbors();
    find_neighbors(X, FORWARD, particles_x_pos, num_particles);
    find_neighbors(X, BACKWARDS, particles_x_pos, num_particles);
    find_neighbors(Y, FORWARD, particles_y_pos, num_particles);
    find_neighbors(Y, BACKWARDS, particles_y_pos, num_particles);
    for (auto particle : particle_system.get_particle_list())
      particle->generate_neighbors_list_by_intersection();

    // Interaction!
    for (auto &particle : particle_system.get_particle_list()) {
      for (auto neighbor : particle->get_neighbors_list()) {
        particle->interact(*neighbor);
      }
    }
    for (auto spring : springs)
      spring->apply_force();

    // Velocity Verlet integration
    particle_system.move_particles(dt, width, height);

    // Data related
    // (TODO: should be moved to ParticleSystem)
    if (step % skip == 0) {
      particle_system.update_trajectory_data(trajectories); // Trajectories
      for (auto particle :
           particle_system.get_particle_list()) { // neighbor matrix
        std::fill(nmat_row.begin(), nmat_row.end(), 0);
        for (auto neighbor : particle->get_neighbors_list())
          nmat_row[neighbor->get_id()] = 1;
        neighbors_matrix.insert(neighbors_matrix.end(), nmat_row.begin(),
                                nmat_row.end());
      }
    }

    // Progress bar update
    progress_perc = (int)((float)step / (float)num_steps * 100);
    bar.set_progress(progress_perc);
    pgtext = "Simulation running: " + std::to_string(step) + "/" +
             std::to_string(num_steps);
    bar.set_option(indicators::option::PostfixText{pgtext});
  }

  // Save data
  std::vector<double> box_size = {width, height};
  // save_data(filename, box_size, num_particles, num_steps, skip,
  trajectories,
      //           masses, radii, neighbors_matrix);
      save_data(filename, box_size, num_particles, num_steps, skip,
                trajectories, masses, radii);

  return 0;
}
