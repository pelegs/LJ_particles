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
  double distances = atof(argv[4]);
  vec2 distances_vec(distances, distances);
  std::string filename = argv[5];

  // randomness!
  srand(time(NULL));
  std::uniform_real_distribution<double> unif_width(0.0, width);
  std::uniform_real_distribution<double> unif_height(0.0, height);
  std::default_random_engine re;
  double rand_x, rand_y;

  ParticleSystem particle_system;
  for (int id=0; id<num_particles; id++) {
    rand_x = unif_width(re);
    rand_y = unif_height(re);
    vec2 pos(rand_x, rand_y);
    particle_system.add_particle(
      new Particle(id, pos, O_, 1.0, 1.0, distances_vec)
    );
  }

  particle_system.calc_new_positions(0.001);
  particle_system.sort_particles_all_directions();
  particle_system.save_data(filename, false, false, true);

  return 0;
}
