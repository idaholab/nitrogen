//* This file is part of nitrogen
//* https://github.com/idaholab/nitrogen
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/nitrogen/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

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
