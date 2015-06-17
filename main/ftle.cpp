#include "ftle.h"
#include "parameter.h"
#include "flowmap.h"
#include "lyapunov.h"
#include "writedata.h"

using namespace std;

int main( int argc, char** argv )
{
cout << "Using " << omp_get_max_threads() << " threads" << endl;

// reading parameter file
auto objpara = make_shared<parameter>(argv[1]);
 
// flowmap calculation
auto objfm = make_shared<flowmap>(objpara);

// ftle calculation from flowmap field
auto objftle = make_shared<lyapunov>(objpara, objfm);

// writing output
auto objwrite = make_shared<writedata>(objpara, objfm, objftle);

return 0;
}