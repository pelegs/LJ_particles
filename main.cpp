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

  // randomness! (note: not even pseudorandom)
  ParticleSystem particle_system(width, height);
  vec2 pos0(width/2.0, height/2.0);
  for (int id = 0; id < num_particles; id++) {
    vec2 pos = glm::sphericalRand(width/2.0);
    vec2 vel = glm::sphericalRand(25.0);
    particle_system.add_particle(
        new Particle(id, pos+pos0, vel, 1.0, 5.0, 50.0));
  }

  for (int step = 0; step < num_steps; step++)
    particle_system.move_particles(dt);

  particle_system.save_data(filename, true, true, true);

  return 0;
}
