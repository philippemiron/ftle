#include "ftle.h"
#include "parameter.h"
#include "flowmap.h"
#include "lyapunov.h"

using namespace std;

int main(int argc, char **argv) {
  cout << "Using " << omp_get_max_threads() << " threads" << endl;

  // reading parameter file
  auto objpara = make_shared<parameter>(argv[1]);

  // flowmap calculation
  auto objfm = make_shared<flowmap>(objpara);
  objfm->calculate_trajectories();

  // ftle calculation from flowmap field
  auto objftle = make_shared<lyapunov>(objpara, objfm);
  objftle->calculate_ftle();
  objftle->output("ftle.dat");
  objftle->outputVTK("ftle.vtk");

  return 0;
}