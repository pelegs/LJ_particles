#include "cnpy.h"
#include "particles.hpp"
#include <set>
#include <string>

template <typename Range, typename Value = typename Range::value_type>
std::string join(Range const);

template <typename container, typename type>
bool in_container(const container &cont, const type &a);

bool particle_in_set(const std::set<Particle *> list, Particle *p);

void save_data(const std::string &filename, const std::vector<double> box_size,
               unsigned long num_particles, unsigned long num_steps,
               unsigned long skip, const std::vector<double> trajectories,
               const std::vector<int> neighbors_matrix,
               const std::vector<double> masses,
               const std::vector<double> radii);
