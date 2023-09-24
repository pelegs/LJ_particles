#include <string>
#include "set"
#include "particles.hpp"

template <typename Range, typename Value = typename Range::value_type>
std::string join(Range const, const char *const);

template <typename container, typename type>
bool in_container(const container &cont, const type &a);
