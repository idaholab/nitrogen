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
// SBTL_call_conv.h
//
#pragma once

#ifdef WIN32
    #ifdef SBTL_EXPORTS
        #define SBTLAPI extern "C" __declspec(dllexport)
    #elif SBTL_IMPORTS
        #define SBTLAPI extern "C" __declspec(dllimport)
    #else
        #define SBTLAPI extern "C"
    #endif
#else
    #define SBTLAPI extern "C"
    #define __stdcall
#endif
