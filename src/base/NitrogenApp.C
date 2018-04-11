#include "NitrogenApp.h"
#include "NitrogenRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "FluidPropertiesApp.h"

// Fluid properties
#include "NitrogenFluidProperties.h"

template <>
InputParameters
validParams<NitrogenApp>()
{
  InputParameters params = validParams<MooseApp>();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

NitrogenApp::NitrogenApp(InputParameters parameters)
  : MooseApp(parameters)
{
  NitrogenApp::registerObjects(_factory);

  FluidPropertiesApp::associateSyntax(_syntax, _action_factory);
}

// External entry point for dynamic application loading
extern "C" void
NitrogenApp__registerApps()
{
  NitrogenApp::registerApps();
}
void
NitrogenApp::registerApps()
{
  registerApp(NitrogenApp);
}

// External entry point for dynamic object registration
extern "C" void
NitrogenApp__registerObjects(Factory & factory)
{
  NitrogenApp::registerObjects(factory);
}

void
NitrogenApp::registerObjects(Factory & factory)
{
  registerUserObject(NitrogenFluidProperties);
}
