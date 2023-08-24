//* This file is part of nitrogen
//* https://github.com/idaholab/nitrogen
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/nitrogen/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseApp.h"

class Factory;

class NitrogenApp : public MooseApp
{
public:
  NitrogenApp(InputParameters parameters);

public:
  static InputParameters validParams();
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};
