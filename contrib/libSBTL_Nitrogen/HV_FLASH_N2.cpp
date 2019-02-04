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
// HV_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_N2.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10

// initial guess from auxiliary splines
extern "C" double __stdcall U_VH_N2_INI_T(double vt, double h);
//
// forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_N2_T(
    double vt, double v, double u, double & p, double & dpdv, double & dpdu, double & dudv);
extern "C" void __stdcall DIFF_P_VU_N2_TT(
    double vt, double u, double & p, double & dpdv, double & dpdu, double & dudv);
//
SBTLAPI int __stdcall FLASH_VH_N2(double v, double h, double & u) throw()
{
  double vt;
  static const double df_h = 1.e-8; //-8   //abs. deviation in h

  double hx, px;
  double dhdu_v;
  double dpdv_u, dpdu_v, dudv_p;

  vt = log(v);

  // calculate initial guess
  u = U_VH_N2_INI_T(vt, h);

  // newtons method
  double f_h = -1.;
  int icount = 0;
  while (fabs(f_h) > df_h)
  {
    DIFF_P_VU_N2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    hx = u + px * v * 1.e3;
    dhdu_v = 1. + dpdu_v * v * 1.e3;
    f_h = hx - h;
    u = u - f_h / dhdu_v;
    if (icount++ > ITMAX)
    {
      u = ERR_VAL;
      return I_ERR;
    }
  }
  return I_OK;
}
//
SBTLAPI int __stdcall FLASH_VH_N2_T(double v, double vt, double h, double & u) throw()
{

  static const double df_h = 1.e-8; //-8   //abs. deviation in h

  double hx, px;
  double dhdu_v;
  double dpdv_u, dpdu_v, dudv_p;

  // calculate initial guess
  u = U_VH_N2_INI_T(vt, h);

  // newtons method
  double f_h = -1.;
  int icount = 0;
  while (fabs(f_h) > df_h)
  {
    DIFF_P_VU_N2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p); // px, transformed derivatives
    hx = u + px * v * 1.e3;
    dhdu_v = 1. + dpdu_v * v * 1.e3;
    f_h = hx - h;
    u = u - f_h / dhdu_v;
    if (icount++ > ITMAX)
    {
      u = ERR_VAL;
      return I_ERR;
    }
  }
  return I_OK;
}
//
SBTLAPI void __stdcall VH_FLASH_DERIV_N2(
    double v, double vt, double u, double & dudv_h, double & dudh_v, double & dvdh_u) throw()
{
  double dpdv_u, dpdu_v, dudv_p;
  double dhdv_u, dhdu_v;
  double p_;

  // derivatives
  DIFF_P_VU_N2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
  dhdv_u = (dpdv_u * v + p_) * 1.e3;
  dhdu_v = 1. + dpdu_v * v * 1.e3;
  dudv_h = -dhdv_u / dhdu_v;
  //
  dudh_v = 1. / dhdu_v;
  dvdh_u = 1. / dhdv_u;
}
//
