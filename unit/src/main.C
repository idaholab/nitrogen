//* This file is part of nitrogen
//* https://github.com/idaholab/nitrogen
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/nitrogen/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"
#include "AppFactory.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "NitrogenApp.h"

PerfLog Moose::perf_log("gtest");

int my_argc;
char ** my_argv;

GTEST_API_ int
main(int argc, char ** argv)
{
  // gtest removes (only) its args from argc and argv - so this  must be before moose init
  testing::InitGoogleTest(&argc, argv);
  my_argc = argc;
  my_argv = argv;

  MooseInit init(argc, argv);
  registerApp(NitrogenApp);
  Moose::_throw_on_error = true;

  return RUN_ALL_TESTS();
}
