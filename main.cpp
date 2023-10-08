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
  double mass = atof(argv[4]);
  double rad = atof(argv[5]);
  double dt = atof(argv[6]);

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
  for (int id = 0; id < num_particles; id++) {
    particle_system.add_particle(new Particle(
        id, grid_points[id], O_, mass, rad, rad * 5.0, sf::Color::Green));
  }
  Particle *p1 = particle_system.get_particle(0);
  Particle *p2 = particle_system.get_particle(1);
  p1->set_color(sf::Color::Blue);
  p2->set_color(sf::Color::White);
  p1->set_vel(50.0, 15.0);

  // Add springs
  double L = glm::distance(p1->get_pos(), p2->get_pos());
  particle_system.add_spring(new Spring(*p1, *p2, 0.1, L));

  // Add walls
  particle_system.add_wall(new Wall(vec2(width, 0.), vec2(0., 0.)));
  particle_system.add_wall(new Wall(vec2(0., 0.), vec2(0., height)));
  particle_system.add_wall(new Wall(vec2(width, 0.), vec2(width, height)));
  particle_system.add_wall(new Wall(vec2(0., height), vec2(width, height)));
  // particle_system.add_wall(
  //     new Wall(vec2(width / 2.0 - 50.0, 100.0), vec2(width / 2.0 + 50.0,
  //     height - 100.0)));

  // Set up SFML window
  int intWidth = (int)width, intHeight = (int)height;
  sf::RenderWindow window(sf::VideoMode(intWidth, intHeight),
                          "SFML graphics test");

  // Setup font
  sf::Font font;
  if (!font.loadFromFile(
          "/usr/share/fonts/truetype/HackNerdFont-Regular.ttf")) {
    std::cout << "Font error" << std::endl;
  }
  sf::Text text;
  int num_chars = 20;
  char buffer[num_chars];
  float fps = 0.0;
  std::snprintf(buffer, num_chars, "fps = ~ %.0f", fps);
  std::string text_str(buffer);
  text.setFont(font); // font is a sf::Font
  text.setString(text_str);
  text.setCharacterSize(100); // in pixels, not points!
  text.setFillColor(sf::Color::Red);
  text.setStyle(sf::Text::Bold);

  // Set up FPS estimation
  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point end;

  // Run
  while (window.isOpen()) {
    sf::Event event;
    start = std::chrono::high_resolution_clock::now();
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed)
        window.close();
    }

    // Dynamics
    particle_system.move_particles(dt, true);

    // ---- Draw ---- //
    window.clear();
    //
    // draw springs
    for (auto spring : particle_system.get_spring_list()) {
      // this entire horrible thing will be replaced with in-spring methods
      vec2 p1_pos = spring->get_particle(0)->get_pos();
      vec2 p2_pos = spring->get_particle(1)->get_pos();
      double p1_rad = spring->get_particle(0)->get_radius();
      double p2_rad = spring->get_particle(1)->get_radius();
      sf::Vertex spring_line[] = {sf::Vertex(sf::Vector2f(p1_pos.x, p1_pos.y)),
                                  sf::Vertex(sf::Vector2f(p2_pos.x, p2_pos.y))};
      window.draw(spring_line, 2, sf::Lines);
    }

    // draw particles
    for (auto particle : particle_system.get_particle_list()) {
      window.draw(particle->get_shape());
    }

    // draw walls
    for (auto wall : particle_system.get_wall_list()) {
      std::array<sf::Vertex, 2> wall_lines =
          wall->get_vertices(); // this is horrible hack that must be corrected
      sf::Vertex line[] = {wall_lines[0], wall_lines[1]};
      window.draw(line, 2, sf::Lines);
    }

    end = std::chrono::high_resolution_clock::now();
    fps = 1.0E9 / (float)std::chrono::duration_cast<std::chrono::nanoseconds>(
                      end - start)
                      .count();
    std::snprintf(buffer, num_chars, "fps = ~ %.0f", fps);
    text_str = std::string(buffer);
    text.setString(text_str);
    window.draw(text);

    window.display();
  }
}
