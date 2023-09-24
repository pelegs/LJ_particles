#ifndef PARTICLES
#define PARTICLES

#include "maths.hpp"
#include <set>

class Particle {
  int id;
  vec2 pos, vel, acc, acc_prev, force;
  double mass, mass_inv, rad;
  std::set<Particle *> neighbors_x, neighbors_y, neighbors;

public:
  Particle();
  Particle(const int &id, const vec2 &pos, const vec2 &vel, const double &mass,
           const double &rad);
  ~Particle();

  // Getters
  int get_id() const;
  vec2 get_pos() const;
  double get_x() const;
  double get_y() const;
  vec2 get_acc_prev() const;
  vec2 get_acc() const;
  vec2 get_force() const;
  double get_mass() const;
  double get_radius() const;
  std::set<Particle *> get_neighbors_list() const;
  std::set<Particle *> get_neighbors_x();
  std::set<Particle *> get_neighbors_y();

  // Setters
  void set_pos(const double &x, const double &y);
  void set_vel(const double &vx, const double &vy);
  void set_mass(const double &m);
  void set_radius(const double &r);

  // Direction between two particles
  vec2 connect(const Particle &p2);
  vec2 look_at(const Particle &p2);

  // Checkers
  bool is_neighbor(Particle *p);
  bool is_neighbor_x(Particle *p);
  bool is_neighbor_y(Particle *p);
  void check_wall_collision(const double &width, const double &height);

  // Force-related stuff
  vec2 LJ_force(const Particle &p2);
  vec2 gravity_force(const Particle &p2);
  void add_force(const vec2 &F);
  void reset_force();
  void interact(const Particle &p2);

  // Velocity Verlet
  void calc_new_pos(const double &dt);
  void calc_acc();
  void calc_new_vel(const double &dt);

  // Neighbors related
  void reset_neighbors();
  void add_neighbor(int axis, Particle *neighbor);
  void generate_neighbors_list_by_intersection();
  std::vector<int> neighbor_ids();
};

#endif // !PARTICLES
