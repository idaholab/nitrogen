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
