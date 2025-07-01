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

Description
    Mean Radiant Temperature (Tmrt) calculation for outdoor thermal comfort
    
    Based on:
    - Thorsson et al. (2007) "Different methods for estimating the mean 
      radiant temperature in an outdoor urban setting"
    - Lindberg et al. (2008) "SOLWEIG 1.0 – Modelling spatial variations 
      of 3D radiant fluxes and mean radiant temperature in complex urban settings"

\*---------------------------------------------------------------------------*/

#ifndef meanRadiantTemperature_H
#define meanRadiantTemperature_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constants
const scalar SIGMA_SB = 5.670374419e-8;  // Stefan-Boltzmann constant [W/m²K⁴]
const scalar EMISSIVITY_HUMAN = 0.97;    // Human body emissivity
const scalar ALBEDO_HUMAN = 0.37;        // Human body albedo (clothing)

// Calculate Tmrt for outdoor environments
inline scalar calculateOutdoorTmrt
(
    const scalar Ta,           // Air temperature [°C]
    const scalar Tsurf,        // Surface temperature [°C]
    const scalar SVF,          // Sky view factor [-]
    const scalar Sdir,         // Direct solar radiation [W/m²]
    const scalar Sdiff,        // Diffuse solar radiation [W/m²]
    const scalar Srefl,        // Reflected solar radiation [W/m²]
    const scalar Ldown,        // Downward longwave radiation [W/m²]
    const scalar Lup           // Upward longwave radiation [W/m²]
)
{
    // Convert temperatures to Kelvin
    scalar Ta_K = Ta + 273.15;
    scalar Tsurf_K = Tsurf + 273.15;
    
    // Sky temperature estimation (simplified)
    // Tsky = Ta * (0.741 + 0.0062 * Tdp)^0.25 where Tdp is dew point
    // For simplification, using clear sky approximation
    scalar Tsky_K = Ta_K * Foam::pow(0.787 + 0.764 * Foam::log((Ta + 273.15)/273.15), 0.25);
    
    // Absorbed shortwave radiation
    scalar Sabs = (1 - ALBEDO_HUMAN) * (Sdir + Sdiff + Srefl);
    
    // Longwave radiation balance
    // For a standing person, projected area factors are considered
    scalar Lnet = EMISSIVITY_HUMAN * (
        SVF * (Ldown - SIGMA_SB * Foam::pow(Tsky_K, 4)) +
        (1 - SVF) * (Lup - SIGMA_SB * Foam::pow(Tsurf_K, 4))
    );
    
    // Total absorbed radiation
    scalar Rabs = Sabs + Lnet;
    
    // Calculate Tmrt
    scalar Tmrt_K = Foam::pow(Rabs / (EMISSIVITY_HUMAN * SIGMA_SB), 0.25);
    
    return Tmrt_K - 273.15;  // Convert back to Celsius
}

// Simplified Tmrt calculation when only basic radiation data available
inline scalar calculateSimplifiedTmrt
(
    const scalar Ta,           // Air temperature [°C]
    const scalar qr,           // Net radiation flux [W/m²]
    const scalar SVF           // Sky view factor [-]
)
{
    if (qr <= 0)
    {
        // If no radiation data, approximate Tmrt ≈ Ta
        return Ta;
    }
    
    // Simple estimation from net radiation
    scalar Ta_K = Ta + 273.15;
    
    // Adjust for sky view factor (outdoor vs indoor)
    scalar adjustmentFactor = 1.0 + 0.3 * SVF;  // Empirical adjustment
    
    scalar Tmrt_K = Foam::pow(
        (qr * adjustmentFactor) / (EMISSIVITY_HUMAN * SIGMA_SB) + 
        Foam::pow(Ta_K, 4), 
        0.25
    );
    
    return Tmrt_K - 273.15;
}

// Calculate Tmrt field for the entire domain
volScalarField calculateMeanRadiantTemperature
(
    const volScalarField& T,
    const volScalarField& qr,
    const volScalarField& SVF,
    const volScalarField* Sdir = nullptr,
    const volScalarField* Sdiff = nullptr
)
{
    volScalarField Tmrt
    (
        IOobject
        (
            "Tmrt_calculated",
            T.time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T.mesh(),
        dimensionedScalar("zero", dimTemperature, 0.0)
    );
    
    forAll(T, cellI)
    {
        // Check if detailed solar radiation data is available
        if (Sdir && Sdiff)
        {
            // Use detailed calculation
            // Note: This would need additional fields for reflected and longwave radiation
            // For now, using simplified version
            Tmrt[cellI] = calculateSimplifiedTmrt(T[cellI], qr[cellI], SVF[cellI]);
        }
        else
        {
            // Use simplified calculation
            Tmrt[cellI] = calculateSimplifiedTmrt(T[cellI], qr[cellI], SVF[cellI]);
        }
    }
    
    return Tmrt;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //