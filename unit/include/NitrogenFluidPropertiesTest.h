#ifndef NITROGENFLUIDPROPERTIESTEST_H
#define NITROGENFLUIDPROPERTIESTEST_H

#include "MooseObjectUnitTest.h"
#include "NitrogenFluidProperties.h"

class NitrogenFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  NitrogenFluidPropertiesTest() : MooseObjectUnitTest("NitrogenApp")
  {
    buildObjects();
  }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("NitrogenFluidProperties");
    _fe_problem->addUserObject("NitrogenFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<NitrogenFluidProperties>("fp");
  }

  const NitrogenFluidProperties * _fp;
};

#endif
