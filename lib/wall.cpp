#include "wall.hpp"

Wall::Wall(const vec2 &p0, const vec2 &p1) {
  this->p0 = p0;
  this->p1 = p1;
  this->dir = p1 - p0;
  this->len = glm::length(this->dir);
  this->dir = glm::normalize(this->dir);
  this->normal = perp2d(this->dir);
}

Wall::Wall(const vec2 &p0, const vec2 &dir, double len) {
  this->p0 = p0;
  this->p1 = p0 + dir * len;
  this->normal = perp2d(this->dir);
}

vec2 Wall::get_p0() const { return this->p0; }
vec2 Wall::get_p1() const { return this->p1; }
vec2 Wall::get_dir() const { return this->dir; }
vec2 Wall::get_normal() const { return this->normal; }
double Wall::get_len() const { return this->len; }
