//* This file is part of nitrogen
//* https://github.com/idaholab/nitrogen
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/nitrogen/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NitrogenApp.h"
#include "NitrogenRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#ifndef SKIP_MODULE_LOAD
#include "ModulesApp.h"
#endif

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

  libmesh_ignore(s);
#ifndef SKIP_MODULE_LOAD
  ModulesApp::registerAllObjects<NitrogenApp>(f, af, s);
#endif
}
