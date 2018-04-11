#include "IAPWS95App.h"
// Moose Includes
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("IAPWS95");

// Begin the main program.
int
main(int argc, char * argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  IAPWS95App::registerApps();

  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("IAPWS95App", argc, argv);

  // Execute the application
  app->run();

  return 0;
}
