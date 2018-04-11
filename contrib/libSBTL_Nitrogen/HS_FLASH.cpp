///////////////////////////////////////////////////////////////////////////
// LibSBTL_vu_N2 - SBTL library for gaseous nitrogen based on:
//
// Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., and Yokozeki, A.:
//
//   "A Reference Equation of State for the Thermodynamic Properties of Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to 2200 MPa,"
//   J. Phys. Chem. Ref. Data, 29(6):1361-1433, 2000.
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
// The Idaho National Laboratory (INL) uses its best efforts to deliver a high-quality software and to verify that the computed information is correct.
// However, INL makes no warranties to that effect, and INL shall not be liable for any damage that may result from errors or omissions in the software.
//
// Version: 0.9.0
//
// HS_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_N2.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10

//initial guess from auxiliary splines
extern "C" void  __stdcall VU_SH_N2_INI(double s, double h, double& vt, double& u);
//
//forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_N2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_N2_T(double vt, double v, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_N2_T(double vt, double v, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_N2_TT(double vt, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_N2_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI int __stdcall HS_FLASH_N2(double h, double s, double& v, double& vt, double& u) throw()
{
    static const double df_h=1.e-8;  //-8   //abs. deviation in h
    static const double df_s=1.e-10; //-8   //abs. deviation in s

    double hx,sx,px,den;
    double dhdv_u, dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;
    double dsdv_u, dsdu_v, dudv_s;

    //calculate initial guesses
    VU_SH_N2_INI(s, h, vt, u);
    v=exp(vt);

    //newtons method
    double f_h=-1.,f_s=-1.;
    int icount=0;
    while(fabs(f_h)>df_h || fabs(f_s)>df_s) {
        DIFF_P_VU_N2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        DIFF_S_VU_N2_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        hx=u+px*v*1.e3;
        dhdv_u=(dpdv_u*v+px*v)*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        f_s=sx-s;
        den=dhdv_u*dsdu_v-dhdu_v*dsdv_u;
        vt=vt+(-dsdu_v*f_h+f_s*dhdu_v)/den;
        u =u +(-f_s*dhdv_u+dsdv_u*f_h)/den;
        v=exp(vt);
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall HS_FLASH_DERIV_N2(double v, double vt, double u, double& dvdh_s, double& dvds_h, double& dhds_v, double& dudh_s, double& duds_h, double& dhds_u) throw()
{
    double dhdv_u, dhdu_v, dudv_h;
    double dpdv_u, dpdu_v, dudv_p;
    double dsdv_u, dsdu_v, dudv_s;
    double p_,s_;

    //derivatives
    DIFF_P_VU_N2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_S_VU_N2_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
    dhdv_u=(dpdv_u*v+p_)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    dudv_h=-dhdv_u/dhdu_v;
    //
    dvdh_s=1./(dhdv_u+dhdu_v*dudv_s);
    dvds_h=1./(dsdv_u+dsdu_v*dudv_h);
    dhds_v=-dvds_h/dvdh_s;
    //
    dudh_s=1./(dhdu_v+dhdv_u/dudv_s);
    duds_h=1./(dsdu_v+dsdv_u/dudv_h);
    dhds_u=-duds_h/dudh_s;
}
//
SBTLAPI void __stdcall HS_PT_FLASH_DERIV_N2(double v, double vt, double u, double& p, double& dpdh_s, double& dpds_h, double& dhds_p, double& t, double& dtdh_s, double& dtds_h, double& dhds_t) throw()
{
    double dhdv_u, dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;
    double dtdv_u, dtdu_v, dudv_t;
    double dsdv_u, dsdu_v, dudv_s;
    double s_;

    //derivatives
    DIFF_P_VU_N2_T(vt, v, u, p, dpdv_u, dpdu_v, dudv_p);
    DIFF_T_VU_N2_T(vt, v, u, t, dtdv_u, dtdu_v, dudv_t);
    DIFF_S_VU_N2_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
    dhdv_u=(dpdv_u*v+p)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    //
    dpdh_s=(dpdv_u*dsdu_v-dpdu_v*dsdv_u)/(dhdv_u*dsdu_v-dhdu_v*dsdv_u);
    dpds_h=(dpdv_u*dhdu_v-dpdu_v*dhdv_u)/(dsdv_u*dhdu_v-dsdu_v*dhdv_u);
    dhds_p=-dpds_h/dpdh_s;
    //
    dtdh_s=(dtdv_u*dsdu_v-dtdu_v*dsdv_u)/(dhdv_u*dsdu_v-dhdu_v*dsdv_u);
    dtds_h=(dtdv_u*dhdu_v-dtdu_v*dhdv_u)/(dsdv_u*dhdu_v-dsdu_v*dhdv_u);
    dhds_t=-dtds_h/dtdh_s;
}
