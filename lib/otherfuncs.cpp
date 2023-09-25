#include <iterator>
#include <iostream>
#include <string>
#include <ostream>
#include <sstream>
#include "otherfuncs.hpp"

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

void save_data(const std::string &filename, const std::vector<double> box_size,
               unsigned long num_particles, unsigned long num_steps,
               unsigned long skip, const std::vector<double> trajectories,
               const std::vector<int> neighbors_matrix,
               const std::vector<double> masses,
               const std::vector<double> radii) {
  cnpy::npz_save(filename, "box_size", &box_size[0], {2}, "w");
  cnpy::npz_save(filename, "neighbors_matrix", &neighbors_matrix[0],
                 {num_steps / skip, num_particles, num_particles}, "a");
  cnpy::npz_save(filename, "trajectories", &trajectories[0],
                 {num_steps / skip, num_particles, 2}, "a");
  cnpy::npz_save(filename, "masses", &masses[0], {num_particles}, "a");
  cnpy::npz_save(filename, "radii", &radii[0], {num_particles}, "a");
  // in each frame i: data[i, :, :].T <-- note the transpose!
}
