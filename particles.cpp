#include "cnpy.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <glm/ext/matrix_transform.hpp>
#include <indicators/progress_bar.hpp>
// #include <indicators/cursor_control.hpp>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <memory>
#include <ostream>
#include <random>
#include <set>
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
const int FORWARD = 1;
const int BACKWARDS = -1;

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

template <typename Range, typename Value = typename Range::value_type>
std::string join(Range const &elements, const char *const delimiter) {
  std::ostringstream os;
  auto b = begin(elements), e = end(elements);

  if (b != e) {
    std::copy(b, prev(e), std::ostream_iterator<Value>(os, delimiter));
    b = prev(e);
  }
  if (b != e) {
    os << *b;
  }

  return os.str();
}

template <typename container, typename type>
bool in_container(const container &cont, const type &a) {
  for (auto vec_element : cont)
    if (vec_element == a)
      return 1;
  return 0;
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
  std::set<Particle *> neighbors_x, neighbors_y, neighbors;

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
  std::set<Particle *> get_neighbors_list() const { return this->neighbors; }
  std::set<Particle *> get_neighbors_x() { return this->neighbors_x; }
  std::set<Particle *> get_neighbors_y() { return this->neighbors_y; }
  bool is_neighbor(Particle *p) { return in_container(this->neighbors, p); }
  bool is_neighbor_x(Particle *p) { return in_container(this->neighbors_x, p); }
  bool is_neighbor_y(Particle *p) { return in_container(this->neighbors_y, p); }

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
      this->neighbors_x.insert(neighbor);
    if (axis == Y)
      this->neighbors_y.insert(neighbor);
  }

  void generate_neighbors_list_by_intersection() {
    this->neighbors.clear();
    set_intersection(this->neighbors_x.begin(), this->neighbors_x.end(),
                     this->neighbors_y.begin(), this->neighbors_y.end(),
                     std::inserter(this->neighbors, this->neighbors.begin()));
    // std::vector<int> neighbors_x_ids = {};
    // for (auto p : this->neighbors_x)
    //   neighbors_x_ids.push_back(p->get_id());
    //
    // std::vector<int> neighbors_y_ids = {};
    // for (auto p : this->neighbors_y)
    //   neighbors_y_ids.push_back(p->get_id());
    //
    // std::vector<int> neighbors_all_ids = {};
    // for (auto p : this->neighbors)
    //   neighbors_all_ids.push_back(p->get_id());
    //
    // std::cerr << "(" << this->id << ") neighbors in x = {"
    //           << join(neighbors_x_ids, ",") << "}, neighbors in y = {"
    //           << join(neighbors_y_ids, ",") << "}, neighbors total = {"
    //           << join(neighbors_all_ids, ",") << "}" << std::endl;
  }

  std::vector<int> neighbor_ids() {
    std::vector<int> ids = {};
    for (auto neighbor : this->neighbors) {
      ids.push_back(neighbor->get_id());
    }
    return ids;
  }

  // void remove_self_from_neighbors_list() {
  //   int index = -1;
  //   for (auto neighbor : this->neighbors) {
  //     ++index;
  //     if (this->id == neighbor->get_id()) {
  //       break;
  //     }
  //   }
  //   this->neighbors.erase(this->neighbors.begin() + index);
  // }

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
               unsigned long num_particles, unsigned long num_steps, unsigned long skip,
               const std::vector<double> trajectories,
               const std::vector<int> neighbors_matrix) {
  cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
  cnpy::npz_save(filename, "neighbors_matrix", &neighbors_matrix[0],
                 {num_steps/skip, num_particles, num_particles}, "a");
  cnpy::npz_save(filename, "trajectories", &trajectories[0],
                 {num_steps/skip, num_particles, 2}, "a");
  // in each frame i: data[i, :, :].T <-- note the transpose!
}

void find_neighbors(const int &axis, const int &dir,
                    const std::vector<Particle *> particle_list,
                    const int &num_particles) {
  for (int id = 0; id < num_particles; id++) {
    for (int i = id + dir; i >= 0 && i < num_particles; i += dir) {
      if (distance1D(particle_list[id]->get_pos()[axis],
                     particle_list[i]->get_pos()[axis]) <= R_CUTOFF)
        particle_list[id]->add_neighbor(axis, particle_list[i]);
      else
        break;
    }
  }
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
  int skip = atoi(argv[5]);
  double dt = atof(argv[6]);
  std::string filename = argv[7];
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
  for (int step = 0; step < num_steps; step++) {
    // Sort x- and y-position vectors
    std::sort(particles_x_pos.begin(), particles_x_pos.end(),
              comapreParticleByXPos);
    std::sort(particles_y_pos.begin(), particles_y_pos.end(),
              comapreParticleByYPos);

    // Create neighbor lists for particles
    for (auto particle : particles)
      particle->reset_neighbors();
    find_neighbors(X, FORWARD, particles_x_pos, num_particles);
    find_neighbors(X, BACKWARDS, particles_x_pos, num_particles);
    find_neighbors(Y, FORWARD, particles_y_pos, num_particles);
    find_neighbors(Y, BACKWARDS, particles_y_pos, num_particles);
    for (auto particle : particles) {
      particle->generate_neighbors_list_by_intersection();
      std::fill(nmat_row.begin(), nmat_row.end(), 0);
      for (auto neighbor : particle->get_neighbors_list())
        nmat_row[neighbor->get_id()] = 1;
      neighbors_matrix.insert(neighbors_matrix.end(), nmat_row.begin(),
                              nmat_row.end());
    }

    // test
    // for (auto particle : particles) {
    //   for (auto neighbor : particle->get_neighbors_list()) {
    //     if (!neighbor->is_neighbor(particle)) {
    //       std::cerr << "Frame: " << step << ", particle " <<
    //       particle->get_id()
    //                 << " is not in neighbor list of its neighbor "
    //                 << neighbor->get_id() << std::endl;
    //       std::cerr << "List of neighbors of particle " << particle->get_id()
    //                 << ": {" << join(particle->neighbor_ids(), ",") << "}"
    //                 << std::endl;
    //       std::cerr << "List of neighbors of particle " << neighbor->get_id()
    //                 << ": {" << join(neighbor->neighbor_ids(), ",") << "}"
    //                 << std::endl;
    //       std::cerr << "-----------------------------------------------"
    //                 << std::endl;
    //     }
    //   }
    // }

    // Interaction!
    for (auto &particle : particles) {
      for (auto neighbor : particle->get_neighbors_list()) {
        particle->interact(*neighbor);
      }
    }

    // Velocity Verlet integration
    move_particles(particles, dt, width, height);

    // Data related
    if (step % skip == 0)
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
  save_data(filename, box_size, num_particles, num_steps, skip, trajectories,
            neighbors_matrix);

  return 0;
}
