#include "cnpy.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <glm/ext/matrix_transform.hpp>
#include <indicators/progress_bar.hpp>
// #include <indicators/cursor_control.hpp>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <ostream>
#include <random>
#include <vector>

// GLM-related
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#define assertm(exp, msg)                                                      \
  assert(((void)msg, exp)) // use (void) to silence unused warnings

/***********************/
/*        Types        */
/***********************/

typedef glm::vec<2, double> vec2;
typedef glm::mat<2, 2, double> mat22;
typedef std::vector<std::vector<int>> int_mat;

/***************************/
/*        Constants        */
/***************************/

// General constants
const double PERCISION = 1.0E-7;
const double inf = std::numeric_limits<double>::infinity();
const double pi = glm::pi<double>();
const double two_pi = 2.0 * pi;
const double half_pi = glm::half_pi<double>();
const double third_pi = pi / 3.0;
const double quarter_pi = half_pi / 2.0;
const double sixth_pi = third_pi / 2.0;
const double sqrt_2 = glm::root_two<double>();
const double one_over_sqrt_2 = 1.0 / sqrt_2;

// Physics-ish constants
const double R_CUTOFF = 30.0;
const double GRAV = 1.0E2;
const double LJ_E = 1.0E6;

// Vector and matrix constants
const vec2 Zero2 = {.0, .0};
const vec2 X_ = {1.0, .0};
const vec2 Y_ = {.0, 1.0};
const vec2 O_ = {.0, .0};
const mat22 I2 = mat22(1.0);

// Row-columns related
const int ROW = 0;
const int COL = 1;
const int X = 0;
const int Y = 1;

/***********************************/
/*        General functions        */
/***********************************/

template <typename T> std::vector<T> arange(T start, T stop, T step = 1) {
  // Equivalent to numpy's arange function
  std::vector<T> values;
  for (T value = start; value < stop; value += step)
    values.push_back(value);
  return values;
}

template <typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in) {
  // Equivalent to numpy's linspace function
  std::vector<double> linspaced;
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);
  if (num == 0) {
    return linspaced;
  }
  if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }
  double delta = (end - start) / (num - 1);
  for (int i = 0; i < num - 1; ++i) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

int *linear_to_square(int id, int num_cols) {
  // Converts 1D index into 2D index.
  // Example: given 5 cells in each row (num_rows), the 1D index 11 is converted
  // into the 2D index [2,1], and the 1D index 37 is converted into the 2D index
  // [7,2].
  // This is the inverse of the function square_to_linear().
  // (NOTE: row counting starts from 0.)
  static int indices[2];
  indices[ROW] = id / num_cols;
  indices[COL] = id % num_cols;
  return indices;
}

int square_to_linear(int i, int j, int num_cols) {
  // Converts 2D index into 1D index.
  // Example: given 5 cells in each row (num_rows), the 2D index [2,1] is
  // converted into the 1D index 11, and the 1D index [7,2] is converted into
  // the 1D index 37. This is the inverse of the function *linear_to_square().
  // (NOTE: row counting starts from 0.)
  return i * num_cols + j;
}

std::vector<int> get_neighboring_indices(int index, int num_rows, int num_cols,
                                         int M, int wrap_x, int wrap_y) {
  // Returns a list of the 1D indices of the neighboring cell with the given 1D
  // index. Neighbors are considered to be all cells that are M={1,2,3,...}
  // cells away from the given cell, either horizontally, vertically or both.
  // The grid can be wrapped both horizontally and vertically.
  int *indices = linear_to_square(index, num_cols);
  int i = indices[0];
  int j = indices[1];
  int I, J;
  std::vector<int> neighbors = {};
  for (int ix = -M; ix <= M; ++ix) {
    I = i + ix;
    if (wrap_x)
      I = ((I % num_rows) + num_rows) % num_rows;
    if (I < 0 || I >= num_rows)
      continue;
    for (int jy = -M; jy <= M; ++jy) {
      J = j + jy;
      if (wrap_y)
        J = ((J % num_cols) + num_cols) % num_cols;
      if (J < 0 || J >= num_cols)
        continue;
      neighbors.push_back(square_to_linear(I, J, num_cols));
    }
  }
  return neighbors;
}

/***********************************/
/*        Physics functions        */
/***********************************/

double distance1D(const double &x, const double &y) { return abs(x - y); }

double U_LJ(double E, double S, double x) {
  // Lennard-Jones potential
  double S_x_6 = std::pow(S / x, 6.0);
  double S_x_12 = std::pow(S_x_6, 2.0);
  return 4.0 * E * (S_x_12 - S_x_6);
}

double F_LJ(double S, double x) {
  // Lennard-Jones force
  double S_6 = std::pow(S, 6.0);
  return -1.0 * 24.0 * LJ_E * S_6 * (std::pow(x, 6.0) - 2 * S_6) /
         std::pow(x, 13.0);
}

/*************************/
/*        Classes        */
/*************************/

class Particle {
  int id;
  int cell_index;
  vec2 pos, vel, acc, acc_prev, force;
  double mass, mass_inv, rad;
  std::vector<Particle *> neighbors_x, neighbors_y, neighbors;

public:
  Particle() {
    id = 0;
    pos = O_;
    vel = O_;
    acc = O_;
    acc_prev = O_;
    force = O_;
    mass = 1.0;
    mass_inv = 1.0;
    rad = 1.0;
    this->reset_neighbors();
    this->reset_cell_index();
  }

  Particle(const int &id, const vec2 &pos, const vec2 &vel, const double &mass,
           const double &rad) {
    this->id = id;
    this->pos = pos;
    this->vel = vel;
    this->mass = mass;
    this->mass_inv = 1.0 / mass;
    this->rad = rad;
    this->reset_neighbors();
    this->reset_cell_index();
  }

  // General getters
  int get_id() const { return this->id; }
  vec2 get_pos() const { return this->pos; }
  double get_x() const { return this->pos[0]; }
  double get_y() const { return this->pos[1]; }
  vec2 get_acc_prev() const { return this->acc_prev; }
  vec2 get_acc() const { return this->acc; }
  vec2 get_force() const { return this->force; }
  double get_mass() const { return this->mass; }
  double get_radius() const { return this->rad; }
  int get_cell_index() const { return this->cell_index; }
  std::vector<Particle *> get_neighbors_list() const { return this->neighbors; }
  std::vector<Particle *> get_neighbors_x() { return this->neighbors_x; }
  std::vector<Particle *> get_neighbors_y() { return this->neighbors_y; }

  // General getters
  void set_pos(const double &x, const double &y) {
    vec2 pos = {x, y};
    this->pos = pos;
  }

  vec2 connect(const Particle &p2) { return p2.get_pos() - this->pos; }
  vec2 look_at(const Particle &p2) { return glm::normalize(this->connect(p2)); }

  // Checking for borders
  void check_wall_collision(const double &width, const double &height) {
    if (pos[X] < this->rad || std::abs(pos[X] - this->rad) > width) {
      this->vel[X] *= -1;
    }
    if (pos[Y] < this->rad || std::abs(pos[Y] - this->rad) > height) {
      this->vel[Y] *= -1;
    }
  }

  // Force-related stuff
  vec2 LJ_force(const Particle &p2) {
    vec2 dir = this->connect(p2);
    double distance = glm::length(dir);
    dir = glm::normalize(dir);
    double F = F_LJ(p2.rad, distance);
    return F * dir;
  }

  vec2 gravity_force(const Particle &p2) {
    vec2 dir = this->connect(p2);
    double distance = glm::length(dir);
    dir = glm::normalize(dir);
    if (distance < p2.get_radius())
      dir *= -1.0;
    double F = GRAV * p2.get_mass() * this->mass / std::pow(distance, 2.0);
    return F * dir;
  }

  void add_force(const vec2 &F) { this->force += F; }
  void reset_force() { this->force = O_; }

  void interact(const Particle &p2) {
    this->add_force(this->LJ_force(p2));
    // this->add_force(this->gravity_force(p2));
  }

  // Velocity Verlet?..
  void calc_new_pos(const double &dt) {
    this->pos += this->vel * dt + 0.5 * this->acc * dt * dt;
  }
  void calc_acc() {
    this->acc_prev = this->acc;
    this->acc = this->force * this->mass_inv;
    this->reset_force();
  }
  void calc_new_vel(const double &dt) {
    vel += 0.5 * (this->acc_prev + this->acc) * dt;
  }

  // Neighbors related
  void reset_neighbors() {
    this->neighbors_x.clear();
    this->neighbors_y.clear();
    this->neighbors.clear();
  }

  void add_neighbor(int axis, Particle *neighbor) {
    if (axis == X)
      this->neighbors_x.push_back(neighbor);
    if (axis == Y)
      this->neighbors_y.push_back(neighbor);
  }

  void set_neighbors_by_intersection() {
    this->neighbors.clear();
    std::set_intersection(this->neighbors_x.begin(), this->neighbors_x.end(),
                     this->neighbors_y.begin(), this->neighbors_y.end(),
                     std::back_inserter(this->neighbors));
  }

  std::vector<int> neighbor_ids() {
    std::vector<int> ids = {};
    for (auto neighbor : this->neighbors) {
      ids.push_back(neighbor->get_id());
    }
    return ids;
  }

  void remove_self_from_neighbors_list() {
    int index = -1;
    for (auto neighbor : this->neighbors) {
      ++index;
      if (this->id == neighbor->get_id()) {
        break;
      }
    }
    this->neighbors.erase(this->neighbors.begin() + index);
  }

  void reset_cell_index() { this->cell_index = -1; }
  void set_cell_index(int index) { this->cell_index = index; }

  // Get particle's info
  void print_data(int print_newline = 0, int print_id = 0,
                  int print_neighbor_ids = 0) const {
    if (print_id)
      std::cout << this->id << ": ";
    std::cout << this->pos[0] << " " << this->pos[1];
    if (print_neighbor_ids) {
      std::cout << "[";
      for (auto neighbor : this->neighbors)
        std::cout << neighbor->get_id() << " ";
    }
    if (print_newline)
      std::cout << std::endl;
    else
      std::cout << " ";
  }
};

/******************************************/
/*        Class-relevant functions        */
/******************************************/

bool comapreParticleByXPos(const Particle *lhs, const Particle *rhs) {
  return lhs->get_x() < rhs->get_x();
}

bool comapreParticleByYPos(const Particle *lhs, const Particle *rhs) {
  return lhs->get_y() < rhs->get_y();
}

void calc_new_positions(std::vector<Particle *> particles, const double &dt) {
  for (auto &p : particles) {
    p->calc_new_pos(dt);
  }
}

void calc_accelerations(std::vector<Particle *> particles) {
  for (auto &p : particles) {
    p->calc_acc();
  }
}

void calc_new_velocities(std::vector<Particle *> particles, const double &dt,
                         const double &width, const double &height) {
  for (auto &p : particles) {
    p->calc_new_vel(dt);
    p->check_wall_collision(width, height);
  }
}

void move_particles(std::vector<Particle *> particles, const double &dt,
                    const double &width, const double &height) {
  calc_new_positions(particles, dt);
  calc_accelerations(particles);
  calc_new_velocities(particles, dt, width, height);
}

void print_particles_data(const std::vector<Particle *> &particles,
                          int print_masses = 1, int print_radii = 1,
                          int print_x = 1, int print_y = 1,
                          int print_cell_index_1D = 1, int num_cols = 0) {
  if (print_masses) {
    std::cout << "# masses:";
    for (auto p : particles)
      std::cout << " " << p->get_mass() << " ";
    std::cout << std::endl;
  }

  if (print_radii) {
    std::cout << "# radii:";
    for (auto p : particles)
      std::cout << " " << p->get_radius() << " ";
    std::cout << std::endl;
  }

  if (print_x) {
    for (auto p : particles)
      std::cout << p->get_x() << " ";
    std::cout << std::endl;
  }

  if (print_y) {
    for (auto p : particles)
      std::cout << p->get_y() << " ";
    std::cout << std::endl;
  }

  if (print_cell_index_1D) {
    for (auto p : particles)
      std::cout << p->get_cell_index() << " ";
    std::cout << std::endl;
  }
}

void append_new_data(const std::vector<Particle *> particles,
                     std::vector<double> &data) {
  for (auto particle : particles) {
    data.push_back(particle->get_x());
    data.push_back(particle->get_y());
  }
}

void save_data(const std::string &filename, const std::vector<double> box_size,
               const std::vector<unsigned long> num_particles,
               unsigned long num_steps, const std::vector<double> trajectories,
               const std::vector<int> neighbors_matrix) {
  cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
  cnpy::npz_save(filename, "num_particles", &num_particles[0], {1}, "a");
  cnpy::npz_save(filename, "neighbors_matrix", &neighbors_matrix[0],
                 {num_steps, num_particles[0], num_particles[0]}, "a");
  cnpy::npz_save(filename, "trajectories", &trajectories[0],
                 {num_steps, num_particles[0], 2}, "a");
  // in each frame i: data[i, :, :].T <-- note the transpose!
}

/**********************/
/*        Main        */
/**********************/

int main(int argc, char *argv[]) {
  // randomness!
  srand(time(NULL));

  // params
  double width = atof(argv[1]);
  double height = atof(argv[2]);
  int num_particles = atoi(argv[3]);
  int num_steps = atoi(argv[4]);
  double dt = atof(argv[5]);
  std::string filename = argv[6];
  filename += ".npz";

  // Neighbor data
  std::vector<int> neighbors_matrix = {};
  std::vector<int> nmat_row(num_particles, 0);

  // init particles
  double x, y;
  std::vector<Particle *> particles;
  for (int i = 0; i < num_particles; i++) {
    vec2 pos(glm::linearRand(0.0, width), glm::linearRand(0.0, height));
    vec2 vel = glm::circularRand(1.0E1);
    particles.push_back(new Particle(i, pos, vel, 1.0, 1.0));
  }

  // Sorted vecs
  std::vector<Particle *> particles_x_pos = {};
  std::vector<Particle *> particles_y_pos = {};
  for (auto particle : particles) {
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

  // Testing neighor finding strategy
  int i, id;
  for (int step = 0; step < num_steps; step++) {
    // Sort x- and y-position vectors
    std::sort(particles_x_pos.begin(), particles_x_pos.end(),
              comapreParticleByXPos);
    std::sort(particles_y_pos.begin(), particles_y_pos.end(),
              comapreParticleByYPos);

    // Create neighbor lists for particle with id (to be generalized soon)
    for (id = 0; id < num_particles; id++) {
      particles[id]->reset_neighbors();
      for (i = id + 1; i < num_particles; i++) {
        if (distance1D(particles[id]->get_x(), particles[i]->get_x()) <=
            R_CUTOFF)
          particles[id]->add_neighbor(X, particles[i]);
        else
          break;
      }
      for (i = id - 1; i > 0; i--) {
        if (distance1D(particles[id]->get_x(), particles[i]->get_x()) <=
            R_CUTOFF)
          particles[id]->add_neighbor(X, particles[i]);
        else
          break;
      }
      for (i = id + 1; i < num_particles; i++) {
        if (distance1D(particles[id]->get_y(), particles[i]->get_y()) <=
            R_CUTOFF)
          particles[id]->add_neighbor(Y, particles[i]);
        else
          break;
      }
      for (i = id - 1; i > 0; i--) {
        if (distance1D(particles[id]->get_y(), particles[i]->get_y()) <=
            R_CUTOFF)
          particles[id]->add_neighbor(Y, particles[i]);
        else
          break;
      }
      particles[id]->set_neighbors_by_intersection();
      std::fill(nmat_row.begin(), nmat_row.end(), 0);
      for (auto neighbor : particles[id]->get_neighbors_list())
        nmat_row[neighbor->get_id()] = 1;
      neighbors_matrix.insert(neighbors_matrix.end(), nmat_row.begin(),
                              nmat_row.end());
    }

    // Velocity Verlet integration
    move_particles(particles, dt, width, height);

    // Data related
    append_new_data(particles, trajectories);

    // Progress bar update
    progress_perc = (int)((float)step / (float)num_steps * 100);
    bar.set_progress(progress_perc);
    pgtext = "Simulation running: " + std::to_string(step) + "/" +
             std::to_string(num_steps);
    bar.set_option(indicators::option::PostfixText{pgtext});
  }

  // Save data
  std::vector<double> box_size = {width, height};
  std::vector<unsigned long> num_particles_vec = {
      static_cast<unsigned long>(num_particles)};
  save_data(filename, box_size, num_particles_vec, num_steps, trajectories,
            neighbors_matrix);

  return 0;
}
