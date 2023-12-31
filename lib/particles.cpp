#include "particles.hpp"
#include "maths.hpp"
#include "otherfuncs.hpp"
#include "physics.hpp"
#include <algorithm>
#include <glm/gtx/closest_point.hpp>
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
  rad_eff = POW_2_1_6;
  bounding_distance = 10.0;
  this->set_bounding_points();
  this->reset_neighbors();

  // Graphics
  this->color = sf::Color::White;
  this->sphere_object = sf::CircleShape(this->rad);
  this->sphere_object.setFillColor(this->color);
}

Particle::Particle(const int &id, const vec2 &pos, const vec2 &vel,
                   const double &mass, const double &rad,
                   const double &bounding_distance, const sf::Color &color, Gradient *colors_gradient) {
  this->id = id;
  this->pos = pos;
  this->vel = vel;
  this->mass = mass;
  this->mass_inv = 1.0 / mass;
  this->rad = rad;
  this->rad_eff = rad * POW_2_1_6;
  this->bounding_distance = bounding_distance;
  this->set_bounding_points();
  this->reset_neighbors();

  // Graphics
  this->color = color;
  this->colors_gradient = colors_gradient;
  this->sphere_object = sf::CircleShape(this->rad * 0.75);
  this->sphere_object.setFillColor(this->color);
}

Particle::~Particle() { this->reset_neighbors(); }

// General getters
int Particle::get_id() const { return this->id; }
vec2 Particle::get_pos() const { return this->pos; }
double Particle::get_x() const { return this->pos[0]; }
double Particle::get_y() const { return this->pos[1]; }
double Particle::get_vx() const { return this->vel[0]; }
double Particle::get_vy() const { return this->vel[1]; }
vec2 Particle::get_acc_prev() const { return this->acc_prev; }
vec2 Particle::get_acc() const { return this->acc; }
vec2 Particle::get_force() const { return this->force; }
double Particle::get_mass() const { return this->mass; }
double Particle::get_radius() const { return this->rad; }
double Particle::get_effective_radius() const { return this->rad_eff; }
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
sf::CircleShape Particle::get_shape() { return this->sphere_object; }

// General setters
void Particle::set_pos(const double &x, const double &y) {
  vec2 pos = {x, y};
  this->pos = pos;
  this->set_bounding_points();
}
void Particle::set_vel(const double &vx, const double &vy) {
  vec2 vel = {vx, vy};
  this->vel = vel;
  this->set_color_from_velocity();
}
void Particle::set_mass(const double &m) {
  this->mass = m;
  this->mass_inv = 1.0 / m;
}
void Particle::set_radius(const double &r) {
  this->rad = r;
  this->rad_eff = r * POW_2_1_6;
}
void Particle::set_bounding_distance(const double &d) {
  this->bounding_distance = d;
  this->set_bounding_points();
}
void Particle::set_bounding_points() {
  vec2 min_point = this->pos - this->bounding_distance;
  vec2 max_point = this->pos + this->bounding_distance;
  this->bounding_points.clear();
  this->bounding_points.push_back(min_point);
  this->bounding_points.push_back(max_point);
}
void Particle::set_color(const sf::Color color) {
  this->color = color;
  this->sphere_object.setFillColor(color);
}

// Direction between two particles
vec2 Particle::connect(const Particle &p2) { return p2.get_pos() - this->pos; }
vec2 Particle::look_at(const vec2 &pt) {
  return glm::normalize(pt - this->pos);
}
vec2 Particle::look_at(const Particle &p2) {
  return glm::normalize(this->connect(p2));
}

// Particle-wall related
vec2 Particle::closest_point_on_wall(const Wall &wall) {
  return glm::closestPointOnLine(this->pos, wall.get_p0(), wall.get_p1());
}
double Particle::distance_to_wall(const Wall &wall) {
  return glm::distance(this->pos, this->closest_point_on_wall(wall)) -
         this->rad;
}
bool Particle::check_collision_with_wall(const Wall &wall, double atol) {
  return this->distance_to_wall(wall) <= atol;
}
void Particle::interact_with_wall(const Wall &wall) {
  // this->vel = glm::reflect(this->vel, wall.get_normal());
  this->add_force(this->WCA_force(wall), MAX_FORCE);
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
  double F = F_LJ(p2.get_radius(), distance);
  return F * dir;
}
vec2 Particle::LJ_force(const Wall &wall) {
  vec2 dir = this->look_at(closest_point_on_wall(wall));
  double distance = glm::length(dir);
  dir = glm::normalize(dir);
  double F = 0.1 * F_LJ(0.75, distance);
  return F * dir;
}
vec2 Particle::WCA_force(const Particle &p2) {
  vec2 dir = this->connect(p2);
  double distance = glm::length(dir);
  dir = glm::normalize(dir);
  double F = F_WCA(p2.get_radius(), distance, this->rad_eff);
  return F * dir;
}
vec2 Particle::WCA_force(const Wall &wall) {
  vec2 dir = this->look_at(closest_point_on_wall(wall));
  double distance = glm::length(dir);
  dir = glm::normalize(dir);
  double F = F_WCA(wall.get_wca_dist(), distance, wall.get_wca_dist_eff());
  return F * dir;
}
vec2 Particle::gravity_force(const Particle &p2) {
  vec2 dir = this->connect(p2);
  double distance2 = glm::length2(dir);
  dir = glm::normalize(dir);
  if (distance2 < p2.get_radius())
    dir *= -1.0;
  double F = GRAV * p2.get_mass() * this->mass / distance2;
  return F * dir;
}
void Particle::add_force(const vec2 &F, double max = -1.0) {
  vec2 actual_F = F;
  if (max >= 0. && glm::length2(F) > std::pow(max, 2.0))
    actual_F = glm::normalize(F) * max;
  this->force += actual_F;
}
void Particle::reset_force() { this->force = O_; }
void Particle::interact_with_particle(const Particle &p2) {
  this->add_force(this->WCA_force(p2), MAX_FORCE);
}

// Velocity Verlet?..
void Particle::calc_new_pos(const double &dt) {
  this->pos += this->vel * dt + 0.5 * this->acc * dt * dt;
  this->set_bounding_points();
  this->update_shape();
}
void Particle::calc_acc() {
  this->acc_prev = this->acc;
  this->acc = this->force * this->mass_inv;
  this->reset_force();
}
void Particle::calc_new_vel(const double &dt) {
  vel += 0.5 * (this->acc_prev + this->acc) * dt;
  this->set_color_from_velocity();
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
  for (auto neighbor : this->neighbors) {
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

// Graphics
void Particle::update_shape() {
  double x = this->pos[X_AX];
  double y = this->pos[Y_AX];
  double r = this->rad;
  this->sphere_object.setPosition(x - rad, y - rad);
}

void Particle::set_color_from_velocity() {
  double vel = glm::length(this->vel);
  this->color = this->colors_gradient->get_color(vel);
  this->sphere_object.setFillColor(this->color);
}

// --------------------------------------------- //

CompareParticlesAABB::CompareParticlesAABB(int ax) { this->axis = ax; }
bool CompareParticlesAABB::operator()(const Particle *p1, const Particle *p2) {
  return p1->get_min_AABB(axis) < p2->get_min_AABB(axis);
}
