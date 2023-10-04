#include "cnpy.h"
#include <SFML/Graphics.hpp>
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
  double r = 7.0;
  for (int id = 0; id < num_particles; id++) {
    particle_system.add_particle(
        new Particle(id, grid_points[id], O_, r, r, 15.0, sf::Color::Green));
  }
  particle_system.get_particle(5)->set_mass(20.0);
  particle_system.get_particle(5)->set_radius(20.0);
  particle_system.get_particle(5)->set_bounding_distance(30.0);
  particle_system.get_particle(5)->set_vel(-140.0, 65.0);

  particle_system.add_wall(new Wall(vec2(width, 0.), vec2(0., 0.)));
  particle_system.add_wall(new Wall(vec2(0., 0.), vec2(0., height)));
  particle_system.add_wall(new Wall(vec2(width, 0.), vec2(width, height)));
  particle_system.add_wall(new Wall(vec2(0., height), vec2(width, height)));
  particle_system.add_wall(
      new Wall(vec2(width / 2.0, 100.0), vec2(width / 2.0, height - 100.0)));

  // Set up SFML window
  sf::Vertex line[] = {sf::Vertex(sf::Vector2f(0.0, 0.0)),
                       sf::Vertex(sf::Vector2f(0.0, 0.0))};
  int intWidth = (int)width, intHeight = (int)height;
  sf::RenderWindow window(sf::VideoMode(intWidth, intHeight),
                          "SFML graphics test");

  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed)
        window.close();
    }

    // Dynamics
    particle_system.move_particles(dt, true);

    // Draw
    window.clear();
    for (auto particle : particle_system.get_particle_list()) {
      window.draw(particle->get_shape());
    }
    for (auto wall : particle_system.get_wall_list()) {
      std::array<sf::Vertex, 2> wall_lines =
          wall->get_vertices(); // this is horrible hack that must be corrected
      sf::Vertex line[] = {wall_lines[0], wall_lines[1]};
      window.draw(line, 2, sf::Lines);
    }
    window.display();
  }
}
