#include <set>
#include <algorithm>
#include "particles.hpp"
#include "maths.hpp"
#include "otherfuncs.hpp"
#include "physics.hpp"

Particle::Particle() {
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

Particle::Particle(const int &id, const vec2 &pos, const vec2 &vel,
                   const double &mass, const double &rad) {
  this->id = id;
  this->pos = pos;
  this->vel = vel;
  this->mass = mass;
  this->mass_inv = 1.0 / mass;
  this->rad = rad;
  this->reset_neighbors();
}

Particle::~Particle() { this->reset_neighbors(); }

// General getters
int Particle::get_id() const { return this->id; }
vec2 Particle::get_pos() const { return this->pos; }
double Particle::get_x() const { return this->pos[0]; }
double Particle::get_y() const { return this->pos[1]; }
vec2 Particle::get_acc_prev() const { return this->acc_prev; }
vec2 Particle::get_acc() const { return this->acc; }
vec2 Particle::get_force() const { return this->force; }
double Particle::get_mass() const { return this->mass; }
double Particle::get_radius() const { return this->rad; }
std::set<Particle *> Particle::get_neighbors_list() const {
  return this->neighbors;
}
std::set<Particle *> Particle::get_neighbors_x() { return this->neighbors_x; }
std::set<Particle *> Particle::get_neighbors_y() { return this->neighbors_y; }

// General setters
void Particle::set_pos(const double &x, const double &y) {
  vec2 pos = {x, y};
  this->pos = pos;
}
void Particle::set_vel(const double &vx, const double &vy) {
  vec2 vel = {vx, vy};
  this->vel = vel;
}
void Particle::set_mass(const double &m) {
  this->mass = m;
  this->mass_inv = 1.0 / m;
}
void Particle::set_radius(const double &r) { this->rad = r; }

// Direction between two particles
vec2 Particle::connect(const Particle &p2) { return p2.get_pos() - this->pos; }
vec2 Particle::look_at(const Particle &p2) {
  return glm::normalize(this->connect(p2));
}

// Checkers
bool Particle::is_neighbor(Particle *p) {
  return in_container(this->neighbors, p);
}
bool Particle::is_neighbor_x(Particle *p) {
  return in_container(this->neighbors_x, p);
}
bool Particle::is_neighbor_y(Particle *p) {
  return in_container(this->neighbors_y, p);
}
void Particle::check_wall_collision(const double &width, const double &height) {
  if (pos[X] < this->rad || std::abs(pos[X] - this->rad) > width) {
    this->vel[X] *= -1;
  }
  if (pos[Y] < this->rad || std::abs(pos[Y] - this->rad) > height) {
    this->vel[Y] *= -1;
  }
}

// Force-related stuff
vec2 Particle::LJ_force(const Particle &p2) {
  vec2 dir = this->connect(p2);
  double distance = glm::length(dir);
  dir = glm::normalize(dir);
  double F = F_LJ(p2.rad, distance);
  return F * dir;
}
vec2 Particle::gravity_force(const Particle &p2) {
  vec2 dir = this->connect(p2);
  double distance = glm::length(dir);
  dir = glm::normalize(dir);
  if (distance < p2.get_radius())
    dir *= -1.0;
  double F = GRAV * p2.get_mass() * this->mass / std::pow(distance, 2.0);
  return F * dir;
}
void Particle::add_force(const vec2 &F) { this->force = this->force + F; }
void Particle::reset_force() { this->force = O_; }
void Particle::interact(const Particle &p2) {
  this->add_force(this->LJ_force(p2));
  // this->add_force(this->gravity_force(p2));
}

// Velocity Verlet?..
void Particle::calc_new_pos(const double &dt) {
  this->pos += this->vel * dt + 0.5 * this->acc * dt * dt;
}
void Particle::calc_acc() {
  this->acc_prev = this->acc;
  this->acc = this->force * this->mass_inv;
  this->reset_force();
}
void Particle::calc_new_vel(const double &dt) {
  vel += 0.5 * (this->acc_prev + this->acc) * dt;
}

// Neighbors related
void Particle::reset_neighbors() {
  this->neighbors_x.clear();
  this->neighbors_y.clear();
  this->neighbors.clear();
}

void Particle::add_neighbor(int axis, Particle *neighbor) {
  if (axis == X)
    this->neighbors_x.insert(neighbor);
  if (axis == Y)
    this->neighbors_y.insert(neighbor);
}

void Particle::generate_neighbors_list_by_intersection() {
  this->neighbors.clear();
  set_intersection(this->neighbors_x.begin(), this->neighbors_x.end(),
                   this->neighbors_y.begin(), this->neighbors_y.end(),
                   std::inserter(this->neighbors, this->neighbors.begin()));
}

std::vector<int> Particle::neighbor_ids() {
  std::vector<int> ids = {};
  for (auto neighbor : this->neighbors) {
    ids.push_back(neighbor->get_id());
  }
  return ids;
}
