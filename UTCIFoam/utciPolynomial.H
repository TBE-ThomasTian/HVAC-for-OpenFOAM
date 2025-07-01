/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    UTCI 6th order polynomial approximation based on:
    Bröde et al. (2012) "Deriving the operational procedure for the 
    Universal Thermal Climate Index (UTCI)"
    International Journal of Biometeorology 56:481-494
    
    The polynomial approximation to the UTCI is defined as:
    UTCI = f(Ta, va, D_Tmrt, Pa)
    where:
    - Ta = air temperature [°C]
    - va = wind speed at 10m height [m/s]
    - D_Tmrt = Tmrt - Ta [°C]
    - Pa = water vapor pressure [kPa]

\*---------------------------------------------------------------------------*/

#ifndef utciPolynomial_H
#define utciPolynomial_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// UTCI polynomial coefficients
// These are example coefficients - the actual values need to be obtained from
// the official UTCI source or publications
const scalar UTCI_OFFSET = 0.0;  // Offset term

// Polynomial structure for UTCI approximation
inline scalar calculateUTCIPolynomial
(
    const scalar Ta,      // Air temperature [°C]
    const scalar va,      // Wind speed at 10m [m/s]
    const scalar D_Tmrt,  // Tmrt - Ta [°C]
    const scalar Pa       // Vapor pressure [kPa]
)
{
    // Convert Pa from hPa to kPa for the polynomial
    scalar Pa_kPa = Pa / 10.0;
    
    // Validate input ranges
    if (Ta < -50 || Ta > 50 || 
        va < 0.5 || va > 30 || 
        D_Tmrt < -30 || D_Tmrt > 70 ||
        Pa_kPa < 0 || Pa_kPa > 5)
    {
        return -9999;  // Invalid input
    }
    
    // UTCI 6th order polynomial approximation
    // Structure: UTCI = Ta + Offset + sum of polynomial terms
    
    scalar utci = Ta;
    
    // First order terms
    utci += 6.07562052e-01;
    utci += -2.27712343e-02 * Ta;
    utci += 8.06470249e-04 * Ta * Ta;
    utci += -1.54271372e-04 * Ta * Ta * Ta;
    utci += -3.24651735e-06 * Ta * Ta * Ta * Ta;
    utci += 7.32602852e-08 * Ta * Ta * Ta * Ta * Ta;
    utci += 1.35959073e-09 * Ta * Ta * Ta * Ta * Ta * Ta;
    
    // Wind speed terms
    utci += -2.25836520e+00 * va;
    utci += 8.80326035e-02 * Ta * va;
    utci += 2.16844454e-03 * Ta * Ta * va;
    utci += -1.53347087e-05 * Ta * Ta * Ta * va;
    utci += -5.72983704e-07 * Ta * Ta * Ta * Ta * va;
    utci += -2.55090145e-09 * Ta * Ta * Ta * Ta * Ta * va;
    
    // Wind speed squared terms
    utci += -7.51269505e-01 * va * va;
    utci += -4.08350271e-03 * Ta * va * va;
    utci += -5.21670675e-05 * Ta * Ta * va * va;
    utci += 1.94544667e-06 * Ta * Ta * Ta * va * va;
    utci += 1.14099531e-08 * Ta * Ta * Ta * Ta * va * va;
    
    // Mean radiant temperature difference terms
    utci += 1.58137256e-01 * D_Tmrt;
    utci += -6.57263143e-05 * Ta * D_Tmrt;
    utci += 2.22697524e-07 * Ta * Ta * D_Tmrt;
    utci += -4.16117031e-08 * Ta * Ta * Ta * D_Tmrt;
    utci += -1.27762753e-09 * Ta * Ta * Ta * Ta * D_Tmrt;
    
    // Vapor pressure terms
    utci += 9.33184380e-01 * Pa_kPa;
    utci += -3.31678953e-03 * Ta * Pa_kPa;
    utci += 1.80309453e-05 * Ta * Ta * Pa_kPa;
    utci += -5.77150105e-08 * Ta * Ta * Ta * Pa_kPa;
    
    // Cross terms (examples - full polynomial has many more)
    utci += 7.22847550e-01 * va * D_Tmrt;
    utci += -8.23033013e-03 * Ta * va * D_Tmrt;
    utci += 1.64473285e-04 * va * va * D_Tmrt;
    utci += -3.85439401e-01 * va * Pa_kPa;
    utci += 3.86783615e-02 * D_Tmrt * Pa_kPa;
    
    // Higher order cross terms would be added here...
    // Note: This is a simplified version. The full UTCI polynomial
    // contains approximately 200 terms up to 6th order
    
    return utci;
}

// Alternative: Offset-based formulation
inline scalar calculateUTCIOffset
(
    const scalar Ta,      // Air temperature [°C]
    const scalar va,      // Wind speed at 10m [m/s]  
    const scalar Tmrt,    // Mean radiant temperature [°C]
    const scalar RH       // Relative humidity [%]
)
{
    // Calculate vapor pressure from RH and Ta
    scalar Es = 6.105 * Foam::exp((17.27 * Ta) / (237.7 + Ta)); // Saturation vapor pressure [hPa]
    scalar Pa = (RH / 100.0) * Es; // Vapor pressure [hPa]
    
    // Calculate D_Tmrt
    scalar D_Tmrt = Tmrt - Ta;
    
    // Call polynomial calculation
    return calculateUTCIPolynomial(Ta, va, D_Tmrt, Pa);
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //