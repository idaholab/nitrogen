#ifndef NITROGENAPP_H
#define NITROGENAPP_H

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
  static void registerObjects(Factory & factory);

protected:
};

#endif /* NITROGENAPP_H */
