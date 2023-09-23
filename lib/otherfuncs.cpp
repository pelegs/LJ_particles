#include <iterator>
#include <iostream>
#include <string>
#include <ostream>
#include <sstream>

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
