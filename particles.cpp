#include "cnpy.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <glm/ext/matrix_transform.hpp>
#include <indicators/progress_bar.hpp>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <memory>
#include <ostream>
#include <random>
#include <set>
#include <vector>

// Own lib
#include "lib/physics.hpp"

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

/*************************/
/*        Classes        */
/*************************/

class Particle {
  int id;
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
  }

  ~Particle() { this->reset_neighbors(); }

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
  std::set<Particle *> get_neighbors_list() const { return this->neighbors; }
  std::set<Particle *> get_neighbors_x() { return this->neighbors_x; }
  std::set<Particle *> get_neighbors_y() { return this->neighbors_y; }
  bool is_neighbor(Particle *p) { return in_container(this->neighbors, p); }
  bool is_neighbor_x(Particle *p) { return in_container(this->neighbors_x, p); }
  bool is_neighbor_y(Particle *p) { return in_container(this->neighbors_y, p); }

  // General setters
  void set_pos(const double &x, const double &y) {
    vec2 pos = {x, y};
    this->pos = pos;
  }
  void set_vel(const double &vx, const double &vy) {
    vec2 vel = {vx, vy};
    this->vel = vel;
  }
  void set_mass(const double &m) {
    this->mass = m;
    this->mass_inv = 1.0 / m;
  }
  void set_radius(const double &r) { this->rad = r; }

  // Direction between two particles
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

  void add_force(const vec2 &F) { this->force = this->force + F; }
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
  }

  std::vector<int> neighbor_ids() {
    std::vector<int> ids = {};
    for (auto neighbor : this->neighbors) {
      ids.push_back(neighbor->get_id());
    }
    return ids;
  }

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

class Spring {
  Particle *p1, *p2;
  double K, L;

public:
  Spring(Particle &p1, Particle &p2, double K, double L) {
    this->p1 = &p1;
    this->p2 = &p2;
    this->K = K;
    if (L == -1.0)
      this->L = glm::distance(p1.get_pos(), p2.get_pos());
    else
      this->L = L;
  }

  void get_data() {
    std::cerr << this->K << " " << this->L << " " << this->p1->get_id() << " "
              << this->p2->get_id() << std::endl;
  }

  double particles_distance() {
    return glm::distance(this->p1->get_pos(), this->p2->get_pos());
  }

  double hook_force(double x) { return this->K * (x - this->L) * -1.0; }

  void apply_force() {
    vec2 dr = this->p2->look_at(*this->p1);
    double x = this->particles_distance();
    vec2 F12 = F_HOOK(this->K, x, this->L) * dr;
    vec2 F21 = -1.0 * F12;
    this->p1->add_force(F12);
    this->p2->add_force(F21);
  }
};

class ParticleSystem {
  int num_particles;
  std::vector<Particle *> particle_list;

public:
  ParticleSystem() {
    this->num_particles = 0;
    this->particle_list.clear();
  }

  ~ParticleSystem() {
    for (auto particle : this->particle_list)
      delete particle;
  }

  // Particle management
  void add_particle(Particle *p) { this->particle_list.push_back(p); }
  void remove_particle() {} // TBW
  Particle *get_particle(int i) { return this->particle_list[i]; }
  std::vector<Particle *> get_particle_list() { return this->particle_list; }

  // Dynamics
  void calc_new_positions(const double &dt) {
    for (auto &p : this->particle_list) {
      p->calc_new_pos(dt);
    }
  }

  void calc_accelerations() {
    for (auto &p : this->particle_list) {
      p->calc_acc();
    }
  }

  void calc_new_velocities(const double &dt, const double &width,
                           const double &height) {
    for (auto &p : this->particle_list) {
      p->calc_new_vel(dt);
      p->check_wall_collision(width, height);
    }
  }

  void move_particles(const double &dt, const double &width,
                      const double &height) {
    calc_new_positions(dt);
    // interaction
    calc_accelerations();
    calc_new_velocities(dt, width, height);
  }

  // Data managment
  void append_new_data(std::vector<double> &data) {
    for (auto particle : this->particle_list) {
      data.push_back(particle->get_x());
      data.push_back(particle->get_y());
    }
  }

  void place_particles_in_grid(const double &width, const double &height,
                               const int &nx, const int &ny) {
    int i = 0;
    int row, col;
    double dx = width / ((double)nx + 1);
    double dy = height / ((double)ny + 1);
    for (auto particle : this->particle_list) {
      row = i / ny + 1;
      col = i % ny + 1;
      particle->set_pos((double)col * dx, (double)row * dy);
      ++i;
    }
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

void save_data(const std::string &filename, const std::vector<double> box_size,
               unsigned long num_particles, unsigned long num_steps,
               unsigned long skip, const std::vector<double> trajectories,
               const std::vector<double> masses,
               const std::vector<double> radii) {
  cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
  // cnpy::npz_save(filename, "neighbors_matrix", &neighbors_matrix[0],
  //                {num_steps / skip, num_particles, num_particles}, "a");
  cnpy::npz_save(filename, "trajectories", &trajectories[0],
                 {num_steps / skip, num_particles, 2}, "a");
  cnpy::npz_save(filename, "masses", &masses[0], {num_particles}, "a");
  cnpy::npz_save(filename, "radii", &radii[0], {num_particles}, "a");
  // in each frame i: data[i, :, :].T <-- note the transpose!
}

// void save_data(const std::string &filename, const std::vector<double>
// box_size,
//                unsigned long num_particles, unsigned long num_steps,
//                unsigned long skip, const std::vector<double> trajectories,
//                const std::vector<double> masses,
//                const std::vector<double> radii) {
//   cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
//   cnpy::npz_save(filename, "trajectories", &trajectories[0],
//                  {num_steps / skip, num_particles, 2}, "a");
//   cnpy::npz_save(filename, "masses", &masses[0], {num_particles}, "a");
//   cnpy::npz_save(filename, "radii", &radii[0], {num_particles}, "a");
//   // in each frame i: data[i, :, :].T <-- note the transpose!
// }

// void save_data(const std::string &filename, const std::vector<double>
// box_size,
//                unsigned long num_particles, unsigned long num_steps,
//                const std::vector<double> trajectories) {
//   cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
//   cnpy::npz_save(filename, "trajectories", &trajectories[0],
//                  {num_steps, num_particles, 2}, "a");
//   // in each frame i: data[i, :, :].T <-- note the transpose!
// }

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
  std::cout << R_CUTOFF << std::endl;
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
  ParticleSystem particle_system;
  for (int i = 0; i < num_particles; i++) {
    particle_system.add_particle(new Particle(i, O_, O_, 1.0, 3.0));
  }
  particle_system.place_particles_in_grid(width, height, 8, 8);
  particle_system.get_particle(20)->set_vel(100.0, 150.0);
  // Set special big particle
  // particle_system.get_particle(49)->set_vel(.0, .0);
  // particle_system.get_particle(49)->set_mass(3.0);
  // particle_system.get_particle(49)->set_radius(5.0);

  // Create mass and radius vectors for saving data
  std::vector<double> masses = {};
  std::vector<double> radii = {};
  for (auto particle : particle_system.get_particle_list()) {
    masses.push_back(particle->get_mass());
    radii.push_back(particle->get_radius());
  }

  // Create springs?
  int bounded_indices[4] = {20, 21, 28, 29};
  std::vector<Spring *> springs = {};
  for (int id1 : bounded_indices)
    for (int id2 : bounded_indices)
      if (id1 != id2) {
        springs.push_back(new Spring(*particle_system.get_particle(id1),
                                     *particle_system.get_particle(id2), 30.0,
                                     -1.0));
      }

  // springs.push_back(new Spring(*particle_system.get_particle(id),
  //                              *particle_system.get_particle(21), 50.0,
  //                              -1.0));

  // Sorted vecs
  std::vector<Particle *> particles_x_pos = {};
  std::vector<Particle *> particles_y_pos = {};
  for (auto particle : particle_system.get_particle_list()) {
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

  // Neighor finding
  for (int step = 0; step < num_steps; step++) {
    // Sort x- and y-position vectors
    std::sort(particles_x_pos.begin(), particles_x_pos.end(),
              comapreParticleByXPos);
    std::sort(particles_y_pos.begin(), particles_y_pos.end(),
              comapreParticleByYPos);

    // Create neighbor lists for particles
    for (auto particle : particle_system.get_particle_list())
      particle->reset_neighbors();
    find_neighbors(X, FORWARD, particles_x_pos, num_particles);
    find_neighbors(X, BACKWARDS, particles_x_pos, num_particles);
    find_neighbors(Y, FORWARD, particles_y_pos, num_particles);
    find_neighbors(Y, BACKWARDS, particles_y_pos, num_particles);
    for (auto particle : particle_system.get_particle_list())
      particle->generate_neighbors_list_by_intersection();

    // Interaction!
    for (auto &particle : particle_system.get_particle_list()) {
      for (auto neighbor : particle->get_neighbors_list()) {
        particle->interact(*neighbor);
      }
    }
    for (auto spring: springs)
      spring->apply_force();

    // Velocity Verlet integration
    particle_system.move_particles(dt, width, height);

    // Data related
    // (TODO: should be moved to ParticleSystem)
    if (step % skip == 0) {
      particle_system.append_new_data(trajectories); // Trajectories
      for (auto particle :
           particle_system.get_particle_list()) { // neighbor matrix
        std::fill(nmat_row.begin(), nmat_row.end(), 0);
        for (auto neighbor : particle->get_neighbors_list())
          nmat_row[neighbor->get_id()] = 1;
        neighbors_matrix.insert(neighbors_matrix.end(), nmat_row.begin(),
                                nmat_row.end());
      }
    }

    // Progress bar update
    progress_perc = (int)((float)step / (float)num_steps * 100);
    bar.set_progress(progress_perc);
    pgtext = "Simulation running: " + std::to_string(step) + "/" +
             std::to_string(num_steps);
    bar.set_option(indicators::option::PostfixText{pgtext});
  }

  // Save data
  std::vector<double> box_size = {width, height};
  // save_data(filename, box_size, num_particles, num_steps, skip,
  trajectories,
      //           masses, radii, neighbors_matrix);
      save_data(filename, box_size, num_particles, num_steps, skip,
                trajectories, masses, radii);

  return 0;
}
