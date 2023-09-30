#ifndef MATHS
#define MATHS

#include <cmath>
#include <glm/ext/matrix_transform.hpp>
#include <vector>

// GLM-related
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#define assertm(exp, msg)                                                      \
  assert(((void)msg, exp)) // use (void) to silence unused warnings


typedef glm::vec<2, double> vec2;
typedef glm::mat<2, 2, double> mat22;
typedef std::vector<std::vector<int>> int_mat;

// General constants
const double PERCISION = 1.0E-7;
const double inf = std::numeric_limits<double>::infinity();
const double pi = glm::pi<double>();
const double two_pi = 2.0 * pi;
const double half_pi = glm::half_pi<double>();
const double third_pi = pi / 3.0;
const double quarter_pi = half_pi / 2.0;
const double sixth_pi = third_pi / 2.0;
const double sqrt_2 = glm::root_two<double>();
const double one_over_sqrt_2 = 1.0 / sqrt_2;

// Vector and matrix constants
const vec2 Zero2 = {.0, .0};
const vec2 X_ = {1.0, .0};
const vec2 Y_ = {.0, 1.0};
const vec2 O_ = {.0, .0};
const mat22 I2 = mat22(1.0);

// Row-columns related
const int ROW = 0;
const int COL = 1;
const int ALL_AXES = -1;
const int X_AX = 0;
const int Y_AX = 1;
const int FORWARD = 1;
const int BACKWARDS = -1;

// Returns a vector of type T elements with values in the interval [start,
// stop), with differnce step
template <typename T> std::vector<T> arange(T, T, T);

// Returns a vector of equaly spaced num_in elements of type T, in the interval
// [start_in, end_in)
template <typename T> std::vector<double> linspace(T, T, int);

// This is used for neighbor finding.
// 1D distance is defined so that there's no need to use expensive squares and
// square root in the case of distances in a single axis.
double distance1D(double, double);

#endif // !MATHS
