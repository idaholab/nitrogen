#pragma once

#include "MooseObjectUnitTest.h"
#include "NitrogenSBTLFluidProperties.h"

class NitrogenSBTLFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  NitrogenSBTLFluidPropertiesTest() : MooseObjectUnitTest("NitrogenApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("NitrogenSBTLFluidProperties");
    _fe_problem->addUserObject("NitrogenSBTLFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<NitrogenSBTLFluidProperties>("fp");
  }

  const NitrogenSBTLFluidProperties * _fp;
};
