#ifndef NITROGENSBTLFLUIDPROPERTIES_H
#define NITROGENSBTLFLUIDPROPERTIES_H

#include "SinglePhaseFluidProperties.h"
#include "NaNInterface.h"

class NitrogenSBTLFluidProperties;

template <>
InputParameters validParams<NitrogenSBTLFluidProperties>();

/**
 * Properties of nitrogen according to Span et al. computed with the SBTL method
 *
 * The following are the units used in libSBTL:
 * - pressure:           MPa
 * - temperature:        K
 * - specific energy:    kJ/kg
 * - density:            kg/m3
 * - speed of sound:     m/s
 * - dynamic viscosity:  Pa-s
 * - conductivity:       W/m-K
 */
class NitrogenSBTLFluidProperties : public SinglePhaseFluidProperties, public NaNInterface
{
public:
  NitrogenSBTLFluidProperties(const InputParameters & parameters);

  virtual Real p_from_v_e(Real v, Real e) const override;
  virtual void p_from_v_e(Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const override;
  virtual Real T_from_v_e(Real v, Real e) const override;
  virtual void T_from_v_e(Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const override;
  virtual Real c_from_v_e(Real v, Real e) const override;
  virtual void c_from_v_e(Real v, Real e, Real & c, Real & dc_dv, Real & dc_de) const override;
  virtual Real cp_from_v_e(Real v, Real e) const override;
  virtual Real cv_from_v_e(Real v, Real e) const override;
  virtual Real mu_from_v_e(Real v, Real e) const override;
  virtual Real k_from_v_e(Real v, Real e) const override;
  virtual Real s_from_v_e(Real v, Real e) const override;
  virtual void s_from_v_e(Real v, Real e, Real & s, Real & ds_dv, Real & ds_de) const override;
  virtual Real s_from_h_p(Real h, Real p) const override;
  virtual void s_from_h_p(Real h, Real p, Real & s, Real & ds_dh, Real & ds_dp) const override;
  virtual Real rho_from_p_s(Real p, Real s) const override;
  virtual void
  rho_from_p_s(Real p, Real s, Real & rho, Real & drho_dp, Real & drho_ds) const override;
  virtual Real e_from_v_h(Real v, Real h) const override;
  virtual void e_from_v_h(Real v, Real h, Real & e, Real & de_dv, Real & de_dh) const override;
  virtual Real rho_from_p_T(Real p, Real T) const override;
  virtual void
  rho_from_p_T(Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const override;
  virtual Real e_from_p_rho(Real p, Real rho) const override;
  virtual void
  e_from_p_rho(Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const override;
  virtual Real e_from_T_v(Real T, Real v) const override;
  virtual void e_from_T_v(Real T, Real v, Real & e, Real & de_dT, Real & de_dv) const override;
  virtual Real p_from_T_v(Real T, Real v) const override;
  virtual void p_from_T_v(Real T, Real v, Real & p, Real & dp_dT, Real & dp_dv) const override;
  virtual Real h_from_T_v(Real T, Real v) const override;
  virtual void h_from_T_v(Real T, Real v, Real & h, Real & dh_dT, Real & dh_dv) const override;
  virtual Real s_from_T_v(Real T, Real v) const override;
  virtual void s_from_T_v(Real T, Real v, Real & s, Real & ds_dT, Real & ds_dv) const override;
  virtual Real cv_from_T_v(Real T, Real v) const override;
  virtual Real h_from_p_T(Real p, Real T) const override;
  virtual void h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const override;
  virtual Real cp_from_p_T(Real p, Real T) const override;
  virtual Real cv_from_p_T(Real p, Real T) const override;
  virtual Real mu_from_p_T(Real p, Real T) const override;
  virtual Real k_from_p_T(Real p, Real T) const override;
  virtual void k_from_p_T(Real p, Real T, Real & k, Real & dk_dp, Real & dk_dT) const override;
  virtual Real p_from_h_s(Real h, Real s) const override;
  virtual void p_from_h_s(Real h, Real s, Real & p, Real & dp_dh, Real & dp_ds) const override;
  virtual Real g_from_v_e(Real v, Real e) const override;
  virtual Real beta_from_p_T(Real p, Real T) const override;
  virtual Real molarMass() const override;

protected:
  /// Conversion factor from Pa to MPa
  const Real _to_MPa;
  /// Conversion factor from MPa to Pa
  const Real _to_Pa;
  /// Conversion factor from J to kJ
  const Real _to_kJ;
  /// Conversion factor from kJ to J
  const Real _to_J;
};

#endif /* NITROGENSBTLFLUIDPROPERTIES_H */
