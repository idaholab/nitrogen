///////////////////////////////////////////////////////////////////////////
// LibSBTL_vu_N2 - SBTL library for gaseous nitrogen based on:
//
// Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., and Yokozeki, A.:
//
//   "A Reference Equation of State for the Thermodynamic Properties of Nitrogen for Temperatures
//   from 63.151 to 1000 K and Pressures to 2200 MPa," J. Phys. Chem. Ref. Data, 29(6):1361-1433,
//   2000.
//
// Lemmon, E.W. and Jacobsen, R.T.:
//
//   "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air"
//   Int. J. Thermophys., 25:21-69, 2004.
//
// Copyright (C) Idaho National Laboratory.
// All rights reserved.
//
// Disclaimer:
// The Idaho National Laboratory (INL) uses its best efforts to deliver a high-quality software and
// to verify that the computed information is correct. However, INL makes no warranties to that
// effect, and INL shall not be liable for any damage that may result from errors or omissions in
// the software.
//
// Version: 0.9.0
//
// PH_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_N2.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10

// initial guess from auxiliary splines
extern "C" void __stdcall VU_HP_N2_INI(double h, double p, double & v, double & u);
//
// forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_N2_T(
    double vt, double v, double u, double & p, double & dpdv, double & dpdu, double & dudv);
extern "C" void __stdcall DIFF_T_VU_N2_T(
    double vt, double v, double u, double & t, double & dtdv, double & dtdu, double & dudv);
extern "C" void __stdcall DIFF_P_VU_N2_TT(
    double vt, double u, double & p, double & dpdv, double & dpdu, double & dudv);
//
SBTLAPI int __stdcall PH_FLASH_N2(double p, double h, double & v, double & vt, double & u) throw()
{
  static const double df_p = 1.e-10; //-8   //rel. deviation in p
  static const double df_h = 1.e-8;  //-8   //abs. deviation in h

  double hx, px, den;
  double dhdv_u, dhdu_v;
  double dpdv_u, dpdu_v, dudv_p;

  // calculate initial guesses
  VU_HP_N2_INI(h, p, v, u);
  vt = log(v);

  // newtons method
  double f_p = -1., f_h = -1., p_inv = 1. / p;
  int icount = 0;
  while (fabs(f_p * p_inv) > df_p || fabs(f_h) > df_h)
  {
    DIFF_P_VU_N2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    hx = u + px * v * 1.e3;
    dhdv_u = (dpdv_u * v + px * v) * 1.e3;
    dhdu_v = 1. + dpdu_v * v * 1.e3;
    f_p = px - p;
    f_h = hx - h;
    den = dhdu_v * dpdv_u - dhdv_u * dpdu_v;
    vt = vt + (-dhdu_v * f_p + f_h * dpdu_v) / den;
    u = u + (-f_h * dpdv_u + dhdv_u * f_p) / den;
    v = exp(vt);
    if (icount++ > ITMAX)
    {
      return I_ERR;
    }
  }
  return I_OK;
}
//
SBTLAPI void __stdcall PH_FLASH_DERIV_N2(double p,
                                         double v,
                                         double vt,
                                         double u,
                                         double & dvdp_h,
                                         double & dvdh_p,
                                         double & dpdh_v,
                                         double & dudp_h,
                                         double & dudh_p,
                                         double & dpdh_u) throw()
{
  double dhdv_u, dhdu_v, dudv_h;
  double dpdv_u, dpdu_v, dudv_p;
  double p_;

  // derivatives
  DIFF_P_VU_N2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
  dhdv_u = (dpdv_u * v + p) * 1.e3;
  dhdu_v = 1. + dpdu_v * v * 1.e3;
  dudv_h = -dhdv_u / dhdu_v;
  //
  dvdp_h = 1. / (dpdv_u + dpdu_v * dudv_h);
  dvdh_p = 1. / (dhdv_u + dhdu_v * dudv_p);
  dpdh_v = -dvdh_p / dvdp_h;
  //
  dudp_h = 1. / (dpdu_v + dpdv_u / dudv_h);
  dudh_p = 1. / (dhdu_v + dhdv_u / dudv_p);
  dpdh_u = -dudh_p / dudp_h;
}
//
SBTLAPI void __stdcall PH_T_FLASH_DERIV_G_N2(double p,
                                             double v,
                                             double vt,
                                             double u,
                                             double & t,
                                             double & dtdp_h,
                                             double & dtdh_p,
                                             double & dpdh_t) throw()
{
  double dhdv_u, dhdu_v;
  double dpdv_u, dpdu_v, dudv_p;
  double dtdv_u, dtdu_v, dudv_t;
  double p_;

  // derivatives
  DIFF_P_VU_N2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
  DIFF_T_VU_N2_T(vt, v, u, t, dtdv_u, dtdu_v, dudv_t);
  dhdv_u = (dpdv_u * v + p) * 1.e3;
  dhdu_v = 1. + dpdu_v * v * 1.e3;
  //
  dtdp_h = (dtdv_u * dhdu_v - dtdu_v * dhdv_u) / (dpdv_u * dhdu_v - dpdu_v * dhdv_u);
  dtdh_p = (dtdv_u * dpdu_v - dtdu_v * dpdv_u) / (dhdv_u * dpdu_v - dhdu_v * dpdv_u);
  dpdh_t = -dtdh_p / dtdp_h;
}
