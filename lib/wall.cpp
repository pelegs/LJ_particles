#include "wall.hpp"
#include <SFML/System/Vector2.hpp>

Wall::Wall(const vec2 &p0, const vec2 &p1) {
  this->p0 = p0;
  this->p1 = p1;
  this->dir = p1 - p0;
  this->len = glm::length(this->dir);
  this->dir = glm::normalize(this->dir);
  this->normal = perp2d(this->dir);
  this->vertices[0] = sf::Vertex(sf::Vector2f(this->p0.x, this->p0.y));
  this->vertices[1] = sf::Vertex(sf::Vector2f(this->p1.x, this->p1.y));
}

Wall::Wall(const vec2 &p0, const vec2 &dir, double len) {
  this->p0 = p0;
  this->p1 = p0 + dir * len;
  this->normal = perp2d(this->dir);
  this->vertices[0] = sf::Vertex(sf::Vector2f(this->p0.x, this->p0.y));
  this->vertices[1] = sf::Vertex(sf::Vector2f(this->p1.x, this->p1.y));
}

vec2 Wall::get_p0() const { return this->p0; }
vec2 Wall::get_p1() const { return this->p1; }
vec2 Wall::get_dir() const { return this->dir; }
vec2 Wall::get_normal() const { return this->normal; }
std::array<sf::Vertex, 2> Wall::get_vertices() const {
  return this->vertices; 
}
