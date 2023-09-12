#include <iostream>
using namespace std;

int main() {
  int b = 3;
  for (int a = -20; a <= 20; a++) {
    int c = ((a % b) + b) % b;
    cout << a << " -> " << c << endl;
  }
  return 0;
}
