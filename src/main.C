//* This file is part of nitrogen
//* https://github.com/idaholab/nitrogen
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/nitrogen/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NitrogenApp.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

PerfLog Moose::perf_log("Nitrogen");

int
main(int argc, char * argv[])
{
  MooseInit init(argc, argv);

  NitrogenApp::registerApps();
  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("NitrogenApp", argc, argv);
  app->run();

  return 0;
}
