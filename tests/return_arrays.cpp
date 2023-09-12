#include <cmath>
#include <ctime>
#include <iostream>

int *return_array(float a, float b){
  static int array[2];
  array[0] = (int)std::floor(a);
  array[1] = (int)std::ceil(b);
  return array;
}

int main() {
  int *arr;
  arr = return_array(10.99999, 6.00001);
  std::cout << *(arr+0) << " " << *(arr+1) << std::endl;
  std::cout << typeid(*arr+0).name() << " " << typeid(*arr+1).name() << std::endl;
  return 0;
}
