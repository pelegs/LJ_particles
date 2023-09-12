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

std::vector<int> neighboring_indices(int index, int num_rows, int num_cols,
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

// /****************************************/
// /*        Classes (headers only)        */
// /****************************************/
//
// class Particle;
// class Cell;
// class Grid;

/*************************/
/*        Classes        */
/*************************/

class Particle {
  int id;
  int cell_index;
  vec2 pos, vel, acc, acc_prev, force;
  double mass, mass_inv, rad;
  std::vector<Particle *> neighbors;

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

  // General getters
  void set_pos(const double &x, const double &y) {
    vec2 pos = {x, y};
    this->pos = pos;
  }

  vec2 connect(const Particle &p2) { return p2.get_pos() - this->pos; }
  vec2 look_at(const Particle &p2) { return glm::normalize(this->connect(p2)); }

  // Checking for borders
  void check_wall_collision(const double &width, const double &height) {
    if (pos[X] < this->rad || std::abs(pos[X] - this->rad) > width)
      this->vel[X] *= -1;
    if (pos[Y] < this->rad || std::abs(pos[Y] - this->rad) > height)
      this->vel[Y] *= -1;
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
  void reset_neighbors() { this->neighbors.clear(); }
  void add_neighbor(Particle *neighbor) { this->neighbors.push_back(neighbor); }
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

class Cell {
  int index;
  std::vector<int *> neighbor_indices;
  std::vector<Particle *> particles;

public:
  void reset() { this->particles.clear(); }

  Cell() {
    this->index = 0;
    this->reset();
  }

  Cell(const int &id) {
    this->index = id;
    this->reset();
  }

  int get_id() const { return this->index; }

  void add_particle(Particle *particle) {
    this->particles.push_back(particle);
    particle->set_cell_index(this->index);
  }

  std::vector<Particle *> get_particles() { return this->particles; }

  void print_particles() {
    for (auto p : this->particles) {
      std::cout << p->get_id() << ": " << glm::to_string(p->get_pos())
                << std::endl;
    }
  }

  int is_empty() const { return this->particles.empty(); }
};

class Grid {
  int num_rows, num_cols;
  double width, height, cell_width, cell_height;
  std::vector<Cell *> cells;

public:
  void init_cells() {
    int cell_index = 0;
    this->cells.clear();
    for (int i = 0; i < this->num_rows; ++i) {
      for (int j = 0; j < this->num_cols; ++j) {
        Cell *new_cell = new Cell(cell_index);
        this->cells.push_back(new_cell);
        cell_index++;
      }
    }
  }

  Grid(const int &num_rows, const int &num_cols, const double &width,
       const double &height) {
    this->num_rows = num_rows;
    this->num_cols = num_cols;
    this->width = width;
    this->height = height;
    this->cell_width = width / (double)num_rows;
    this->cell_height = height / (double)num_cols;
    this->init_cells();
  }

  // Getters
  int get_num_rows() { return this->num_rows; }
  int get_num_cols() { return this->num_cols; }
  double get_width() { return this->width; }
  double get_height() { return this->height; }
  double get_cell_width() { return this->cell_width; }
  double get_cell_height() { return this->cell_height; }

  void print_cell_info() {
    for (auto cell : this->cells) {
      // print using cell function
    }
  }

  void reset() {
    for (auto cell : this->cells) {
      cell->reset();
    }
  }

  int get_index_from_pos(const vec2 &pos) const {
    int cell_index, row_index, col_index;
    if (pos[0] < 0.0 || pos[1] < 0.0 || pos[0] >= this->width ||
        pos[1] >= this->height) {
      return -1;
    }
    row_index = (int)std::floor((pos[Y] / this->height) * (this->num_rows));
    col_index = (int)std::floor((pos[X] / this->width) * (this->num_cols));

    int index1D = square_to_linear(row_index, col_index, this->num_cols);
    // std::cout << glm::to_string(pos) << " -> ";
    // std::cout << row_index << ", " << col_index << " -> " << index1D <<
    // std::endl;
    return index1D;
  }

  void add_particles(std::vector<Particle *> particles) {
    int cell_index;
    for (auto particle : particles) {
      cell_index = this->get_index_from_pos(particle->get_pos());
      // std::cout << glm::to_string(particle->get_pos()) << " -> " <<
      // cell_index
      //           << std::endl;
      this->cells[cell_index]->add_particle(particle);
    }
  }

  std::vector<Particle *> get_particles_from_cell(int index) {
    return this->cells[index]->get_particles();
  }

  std::vector<Particle *> get_particles_from_neighboring_cels(int index, int M,
                                                              int wrap_x,
                                                              int wrap_y) {
    std::vector<Particle *> neighbors = {};
    std::vector<Particle *> particles_in_cell = {};
    std::vector<int> neighbor_cell_1D_indices = neighboring_indices(
        index, this->num_rows, this->num_cols, M, wrap_x, wrap_y);
    for (auto idx : neighbor_cell_1D_indices) {
      particles_in_cell = this->cells[idx]->get_particles();
      for (auto particle : particles_in_cell) {
        neighbors.push_back(particle);
      }
    }
    return neighbors;
  }

  void generate_neighbor_list(Particle *particle, int M, int wrap_x,
                              int wrap_y) {
    particle->reset_neighbors();
    std::vector<Particle *> neighbors =
        this->get_particles_from_neighboring_cels(particle->get_cell_index(), M,
                                                  wrap_x, wrap_y);
    for (auto neighbor : neighbors) {
      if (neighbor->get_id() != particle->get_id())
        particle->add_neighbor(neighbor);
    }
  }

  std::vector<vec2> lattice_points() {
    double dx = this->width / (double)(this->num_cols + 1);
    double dy = this->height / (double)(this->num_rows + 1);
    std::vector<vec2> points = {};
    for (int i = 1; i <= this->num_cols; ++i) {
      for (int j = 1; j <= this->num_rows; ++j) {
        points.push_back(vec2(i * dx, j * dy));
      }
    }
    return points;
  }
};

/******************************************/
/*        Class-relevant functions        */
/******************************************/

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
               unsigned long num_steps,
               const std::vector<double> data) {
  cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
  cnpy::npz_save(filename, "num_particles", &num_particles[0], {1}, "a");
  cnpy::npz_save(filename, "trajectories", &data[0],
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

  // test
  int approx_sqrt_num_particles = (int)ceil(std::sqrt(num_particles));
  Grid grid(approx_sqrt_num_particles, approx_sqrt_num_particles, width,
            height);
  std::vector<vec2> points = grid.lattice_points();

  // init particles
  double x, y;
  int index1D;
  std::vector<Particle *> particles;
  for (int i = 0; i < num_particles; i++) {
    vec2 vel = glm::circularRand(1.0E0);
    particles.push_back(new Particle(i, points[i], vel, 1.0, 1.0E-1));
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
  std::vector<double> data = {};
  double time = 0.0;
  for (int step = 0; step < num_steps; step++) {
    for (auto &p1 : particles) {
      for (auto p2 : particles) {
        if (p1->get_id() != p2->get_id()) {
          p1->interact(*p2);
        }
      }
    }

    for (auto &p : particles) {
      p->calc_new_pos(dt);
    }
    for (auto &p : particles) {
      p->calc_acc();
    }
    for (auto &p : particles) {
      p->calc_new_vel(dt);
      p->check_wall_collision(width, height);
    }

    append_new_data(particles, data);
    time += dt;

    progress_perc = (int)((float)step / (float)num_steps * 100);
    bar.set_progress(progress_perc);
    pgtext = "Simulation running: " + std::to_string(step) + "/" + std::to_string(num_steps);
    bar.set_option(indicators::option::PostfixText{pgtext});
  }

  // Save data
  std::vector<double> box_size = {width, height};
  std::vector<unsigned long> num_particles_vec = {
      static_cast<unsigned long>(num_particles)};
  save_data(filename, box_size, num_particles_vec,
            num_steps, data);

  return 0;
}
