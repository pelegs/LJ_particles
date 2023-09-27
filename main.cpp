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

/**********************/
/*        Main        */
/**********************/

int main(int argc, char *argv[]) {
  // params
  double width = atof(argv[1]);
  double height = atof(argv[2]);
  int num_particles = atoi(argv[3]);
  double max_BB_size = atof(argv[4]);
  double max_radius = atof(argv[5]);
  std::string filename = argv[6];

  // randomness! (note: not even pseudorandom)
  std::random_device r;
  std::uniform_real_distribution<double> unif_width(0.0, width);
  std::uniform_real_distribution<double> unif_height(0.0, height);
  std::uniform_real_distribution<double> unif_radii(0.5, max_radius);
  std::uniform_real_distribution<double> unif_bd(1.0, max_BB_size);
  std::default_random_engine re;
  double rand_x, rand_y, rand_rad, rand_bd;

  ParticleSystem particle_system;
  for (int id=0; id<num_particles; id++) {
    rand_x = unif_width(re);
    rand_y = unif_height(re);
    rand_rad = unif_radii(re);
    rand_bd = unif_bd(re);
    vec2 pos(rand_x, rand_y);
    particle_system.add_particle(
      new Particle(id, pos, O_, 1.0, rand_rad, rand_rad+rand_bd)
    );
  }

  particle_system.calc_new_positions(0.001);
  particle_system.reset_neighbors();
  particle_system.sort_particles_all_directions();
  particle_system.assign_neighbors();
  particle_system.update_neighbors_matrix();
  particle_system.save_data(filename, true, true, true);

  return 0;
}
