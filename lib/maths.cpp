#include <cmath>
#include <vector>

template <typename T> std::vector<T> arange(T start, T stop, T step = 1) {
  // Equivalent to numpy's arange function
  std::vector<T> values;
  for (T value = start; value < stop; value += step)
    values.push_back(value);
  return values;
}

template <typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in) {
  // Equivalent to numpy's linspace function
  std::vector<double> linspaced;
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);
  if (num == 0) {
    return linspaced;
  }
  if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }
  double delta = (end - start) / (num - 1);
  for (int i = 0; i < num - 1; ++i) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

double distance1D(double x, double y) {
  return std::abs(x - y);
}
