#include "particles.hpp"
#include "maths.hpp"
#include "otherfuncs.hpp"
#include "physics.hpp"
#include <algorithm>
#include <glm/gtx/string_cast.hpp>
#include <set>

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
  bounding_distance = 10.0;
  this->set_bounding_points();
  this->reset_neighbors();
}

Particle::Particle(const int &id, const vec2 &pos, const vec2 &vel,
                   const double &mass, const double &rad,
                   const double &bounding_distance) {
  this->id = id;
  this->pos = pos;
  this->vel = vel;
  this->mass = mass;
  this->mass_inv = 1.0 / mass;
  this->rad = rad;
  this->bounding_distance = bounding_distance;
  this->set_bounding_points();
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
double Particle::get_bounding_distance() const {
  return this->bounding_distance;
}
double Particle::get_min_AABB(int axis) const {
  return this->bounding_points[MIN_BB][axis];
}
double Particle::get_max_AABB(int axis) const {
  return this->bounding_points[MAX_BB][axis];
}
std::set<Particle *> Particle::get_neighbors_list() const {
  return this->neighbors;
}
std::set<Particle *> Particle::get_neighbors_x() { return this->neighbors_x; }
std::set<Particle *> Particle::get_neighbors_y() { return this->neighbors_y; }

// General setters
void Particle::set_pos(const double &x, const double &y) {
  vec2 pos = {x, y};
  this->pos = pos;
  this->set_bounding_points();
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
void Particle::set_bounding_points() {
  vec2 min_point = this->pos - this->bounding_distance;
  vec2 max_point = this->pos + this->bounding_distance;
  this->bounding_points.clear();
  this->bounding_points.push_back(min_point);
  this->bounding_points.push_back(max_point);
}

// Direction between two particles
vec2 Particle::connect(const Particle &p2) { return p2.get_pos() - this->pos; }
vec2 Particle::look_at(const Particle &p2) {
  return glm::normalize(this->connect(p2));
}

// Checkers
void Particle::check_wall_collision(const double &width, const double &height) {
  if (pos[X_AX] < this->rad || std::abs(pos[X_AX] - this->rad) > width) {
    this->vel[X_AX] *= -1;
  }
  if (pos[Y_AX] < this->rad || std::abs(pos[Y_AX] - this->rad) > height) {
    this->vel[Y_AX] *= -1;
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
  this->set_bounding_points();
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
  this->neighbors.clear();
  this->neighbors_x.clear();
  this->neighbors_y.clear();
}

void Particle::add_neighbor(int axis, Particle *neighbor) {
  if (axis == ALL_AXES)
    this->neighbors.insert(neighbor);
  else if (axis == X_AX)
    this->neighbors_x.insert(neighbor);
  else if (axis == Y_AX)
    this->neighbors_y.insert(neighbor);
}

void Particle::generate_neighbors_list_by_intersection() {
  set_intersection(this->neighbors_x.begin(), this->neighbors_x.end(),
                   this->neighbors_y.begin(), this->neighbors_y.end(),
                   std::inserter(this->neighbors, this->neighbors.begin()));
  for (auto neighbor:this->neighbors) {
    neighbor->add_neighbor(ALL_AXES, this);
  }
}

std::vector<int> Particle::neighbor_ids() {
  std::vector<int> ids = {};
  for (auto neighbor : this->neighbors) {
    ids.push_back(neighbor->get_id());
  }
  return ids;
}

CompareParticlesAABB::CompareParticlesAABB(int ax) { this->axis = ax; }
bool CompareParticlesAABB::operator()(const Particle *p1, const Particle *p2) {
  return p1->get_min_AABB(axis) < p2->get_min_AABB(axis);
}
