!template load file=sqa/app_sdd.md.template app=Nitrogen category=nitrogen

Nitrogen properties are computed within a fluid property user object. Fluid properties user objects derive
from base classes defined in the MOOSE fluid properties module. These base classes define interfaces that other
objects, such as Materials, can leverage to make the fluid property user object retrieve the desired fluid properties,
e.g. density or specific heat.

Nitrogen properties do not include phase changes, so the derived class naturally
derives from the [SinglePhaseFluidProperties.md] base class in the Fluid Properties module.
As such most interfaces that are defined have a similar signature to the one below. Alternative signatures are available for using automatic differentiation or using the specific volume and specific enthalpy variables (also referred to as the conservative variables set).

```
rho_from_p_T(Real pressure, Real temperature)
```

As with any fluid property, the ranges of validity of these fluid properties used are finite. They are documented in the class documentation
for the fluid property user object for SBTL Nitrogen fluid properties. {{app}} may error if these properties are queried outside
of their range of validity, depending on the user-selected out-of-bounds behavior.
