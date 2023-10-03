#include "maths.hpp"
#include <vector>

class Wall {
  vec2 p0, p1, dir, normal;
  double len;

public:
  Wall(const vec2 &p0, const vec2 &p1);
  Wall(const vec2 &p0, const vec2 &dir, double len);

  // getters
  vec2 get_p0() const;
  vec2 get_p1() const;
  vec2 get_dir() const;
  vec2 get_normal() const;
  double get_len() const;
};
