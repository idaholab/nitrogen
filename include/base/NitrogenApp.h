#pragma once

#include "MooseApp.h"

class Factory;
class NitrogenApp;

template <>
InputParameters validParams<NitrogenApp>();

class NitrogenApp : public MooseApp
{
public:
  NitrogenApp(InputParameters parameters);

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

protected:
};
