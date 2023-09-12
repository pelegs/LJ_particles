#include <iostream>
#include <random>

int *linear_to_square(int id, int nx) {
  static int indices[2];
  indices[0] = id / nx;
  indices[1] = id % nx;
  return indices;
}

int square_to_linear(int i, int j, int nx) { return i * nx + j; }

int main(int argc, char *argv[]) {
  const int id_low = 0;
  const int id_high = 100;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> id_distrib(id_low, id_high);
  int id = id_distrib(gen);

  const int nxy_low = 4;
  const int nxy_high = 20;
  std::uniform_int_distribution<> nxy_distrib(nxy_low, nxy_high);
  int nx, ny;
  int found = 0;
  while (!found) {
    nx = nxy_distrib(gen);
    ny = nxy_distrib(gen);
    if (nx * ny > id)
      found = true;
  }

  int *indices = linear_to_square(id, nx);
  int id_back = square_to_linear(indices[0], indices[1], nx);
  std::cout << "nx=" << nx << ", ny=" << ny << std::endl;
  std::cout << id << " -> [" << *(indices + 0) << ", " << *(indices + 1)
            << "] -> " << id_back << std::endl;

  return 0;
}
