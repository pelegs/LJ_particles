#include "gradients.hpp"
#include <iostream>

Gradient::Gradient(const sf::Color &start_color, const sf::Color &end_color,
                   const int &num_steps, const double &start_value = 0.0,
                   const double &end_value = 1.0) {
  this->start_color = start_color;
  this->end_color = end_color;
  this->num_steps = num_steps;
  this->start_value = start_value;
  this->end_value = end_value;

  double dR = (int)((end_color.r - start_color.r) / (double)num_steps);
  double dG = (int)((end_color.g - start_color.g) / (double)num_steps);
  double dB = (int)((end_color.b - start_color.b) / (double)num_steps);
  sf::Color color;
  for (int step = 0; step < num_steps; ++step) {
    color.r += dR;
    color.g += dG;
    color.b += dB;
    this->colors.push_back(color);
  }
}

sf::Color Gradient::get_color(const double &value = 0.0) const {
  if (value <= this->start_value)
    return this->start_color;
  if (value >= this->end_value)
    return this->end_color;
  int index =
      (int)((this->end_value - this->start_value) * value / this->num_steps);
  return this->colors[index];
}
