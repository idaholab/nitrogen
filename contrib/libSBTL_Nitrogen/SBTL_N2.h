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
// SBTL_N2.h
//
#pragma once
//
//-----------------------------------------------------------------------------
// return values (error flags)
//-----------------------------------------------------------------------------
//
#define ERR_VAL -1.e3
#define I_OK  0
#define I_ERR 1

//-----------------------------------------------------------------------------
// struct states
//-----------------------------------------------------------------------------
//
                    //known properties:
#define STR_ERR 0   // - none (needs to be initialized with ireg-call) - currently LibSBTL_vu_N2 considered gaseous nitrogen only (single phase)
#define STR_PDP 1   // - region
                    // - auxiliary variables for property calculations
                    // - auxiliary variables for derivatives in the single-phase region
//
//-----------------------------------------------------------------------------
// structs to handle preliminary results and auxiliary variables:
//-----------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////
// - struct to be used by ireg_vu_SBTL_N2
// - this is also used by ireg_pv_SBTL_N2 and ireg_ps_SBTL_N2
///////////////////////////////////////////////////////////////////////////////
//
typedef struct _STR_vu_SBTL_N2 {
//
    double v_;      //values below computed for this
    double u_;      //state point
    double p_;      //check p_,v_ for given (p,v), p_,s_ for given (p,s), and p_,h_ for given (p,h)
    double h_;      //check p_,h_ for given (p,h) and h_,s_ for given (h,s)
    double s_;      //check p_,s_ for given (p,s) and h_,s_ for given (h,s)
//
    double vt;      //transformed volume        (gas phase only)
// constructor
    _STR_vu_SBTL_N2() { reset();}
// reset
    void reset() {
        v_      =ERR_VAL;
        u_      =ERR_VAL;
        p_      =ERR_VAL;
        h_      =ERR_VAL;
        s_      =ERR_VAL;
        //
        vt      =ERR_VAL;
    }
//
    int GetState(double v, double u) {
        if(v_==v && u_==u) {
            return STR_PDP;
        } else {
            return STR_ERR;
        }
    }
//
    int GetStatePV(double p, double v) {
        if(p_==p && v_==v) {
            return STR_PDP;
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStatePS(double p, double s) {
        if(p_==p && s_==s) {
            return STR_PDP;
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStatePH(double p, double h) {
        if(p_==p && h_==h) {
            return STR_PDP;
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStateHS(double h, double s) {
        if(h_==h && s_==s) {
            return STR_PDP;
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStateVH(double v, double h) {
        if(v_==v && h_==h) {
            return STR_PDP;
        } else {
            reset();
            return STR_ERR;
        }
    }
} STR_vu_SBTL_N2;
