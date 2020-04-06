#include "NitrogenApp.h"
#include "NitrogenRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "ModulesApp.h"

InputParameters
NitrogenApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

registerKnownLabel("NitrogenApp");

NitrogenApp::NitrogenApp(InputParameters parameters) : MooseApp(parameters)
{
  NitrogenApp::registerAll(_factory, _action_factory, _syntax);
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

// External entry point for dynamic application loading
extern "C" void
NitrogenApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  NitrogenApp::registerAll(f, af, s);
}

void
NitrogenApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"NitrogenApp"});
  Registry::registerActionsTo(af, {"NitrogenApp"});

  ModulesApp::registerAll(f, af, s);
}
