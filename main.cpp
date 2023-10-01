#include "cnpy.h"
#include <cassert>
#include <glm/ext/matrix_transform.hpp>
#include <indicators/progress_bar.hpp>
#include <iterator>
#include <random>

// Own lib
#include "lib/maths.hpp"
#include "lib/otherfuncs.hpp"
#include "lib/particle_system.hpp"
#include "lib/particles.hpp"
#include "lib/physics.hpp"
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
  int num_steps = atoi(argv[4]);
  double dt = atof(argv[5]);
  double max_BB_size = atof(argv[6]);
  double max_radius = atof(argv[7]);
  std::string filename = argv[8];

  // init randomness?
  srand(time(NULL));

  ParticleSystem particle_system(width, height);
  double x, y;
  vec2 v = O_;
  for (int id = 0; id < num_particles; id++) {
    x = ((id % 5)+1) * width / 6;
    y = ((id / 5)+1) * height / 6;
    v = glm::diskRand(250.0);
    particle_system.add_particle(
        new Particle(id, vec2(x, y), v, 1.0, 10.0, 20.0));
  }

  for (int step = 0; step < num_steps; step++)
    particle_system.move_particles(dt, true);

  particle_system.save_data(filename, true, true, true, true, false);

  return 0;
}
