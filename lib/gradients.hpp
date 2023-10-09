#include <SFML/Graphics.hpp>

class Gradient {
  sf::Color start_color, end_color;
  double start_value, end_value;
  int num_steps;
  std::vector<sf::Color> colors;

public:
  Gradient(const sf::Color &start_color, const sf::Color &end_color,
           const int &num_steps, const double &start_value,
           const double &end_value);
  sf::Color get_color(const double &value) const;
};
