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
  int root_num_particles = atoi(argv[3]);
  int num_particles = root_num_particles * root_num_particles;
  int num_steps = atoi(argv[4]);
  double dt = atof(argv[5]);
  double max_radius = atof(argv[6]);
  std::string filename = argv[7];

  // init randomness?
  srand(time(NULL));

  // Grid
  std::vector<vec2> grid_points = {};
  int row, col;
  for (int i = 0; i < num_particles; i++) {
    col = (i % root_num_particles + 1) * width / (root_num_particles + 1);
    row = (i / root_num_particles + 1) * height / (root_num_particles + 1);
    grid_points.push_back(vec2(col, row));
  }

  ParticleSystem particle_system(width, height);
  double x, y;
  double r=2.0;
  for (int id = 0; id < num_particles; id++) {
    particle_system.add_particle(
        new Particle(id, grid_points[id], O_, r*5.0, r, r+3.0));
  }
  particle_system.get_particle(25)->set_mass(100.0);
  particle_system.get_particle(25)->set_radius(10.0);
  particle_system.get_particle(25)->set_bounding_distance(17.0);
  particle_system.get_particle(25)->set_vel(40.0, -65.0);

  for (int step = 0; step < num_steps; step++)
    particle_system.move_particles(dt, true);

  particle_system.save_data(filename, true, true, false, false, false);

  return 0;
}
