#include <signal.h>
#include "pch.h"
#include "staticlib.h"


////////////////////////////////////////////////////////////////////////////



void signal_callback(int sig) {
  fmt::print(
      fmt::fg(fmt::color::red), "Interrupt signal [{0}] received, terminating program.\n\n\n\n\n\n\n\n", sig);
  exit(sig);
}

int main() {
  signal(SIGINT, signal_callback);


  fmt::print("Pardiso Works?:{0:}\n", test_pardiso());

  return 0;
}


PYBIND11_MODULE(_asset_skel, m) {


  signal(SIGINT, signal_callback);
  m.def("main", &main);


}
