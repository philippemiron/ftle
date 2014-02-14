#include "ftle.h"
#include "parameter.h"
#include "flowmap.h"
#include "lyapunov.h"
#include "writedata.h"

using namespace std;

int main( int argc, char** argv )
{

// output the number of thread created by openMP
cout << "Using " << omp_get_max_threads() << " threads" << endl;

// reading parameter file
const char* fichierconfig=argv[1];
parameter objpara(fichierconfig);
 
// flowmap calculation
flowmap objfm(&objpara);

// ftle calculation from flowmap field 
lyapunov objftle(&objpara, &objfm);

// writing output
writedata objout(&objpara, &objfm, &objftle);

return 0;

}
