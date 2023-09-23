#include <string>

template <typename Range, typename Value = typename Range::value_type>
std::string join(Range const, const char *const);

template <typename container, typename type>
bool in_container(const container, const type);
