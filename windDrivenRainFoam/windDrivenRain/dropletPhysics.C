/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dropletPhysics.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::volScalarField Foam::dropletPhysics::evaporationRate
(
    const volScalarField& diameter,
    const volScalarField& temperature,
    const volScalarField& relativeHumidity,
    const volScalarField& pressure,
    const volVectorField& Urel
)
{
    // Ranz-Marshall correlation for mass transfer
    volScalarField Re = mag(Urel)*diameter*1.2/1.8e-5; // Approximate air properties
    volScalarField Sc(0.6); // Schmidt number for water vapor in air
    volScalarField Sh = 2.0 + 0.6*sqrt(Re)*pow(Sc, 0.333);
    
    // Vapor pressure difference
    volScalarField psat = 610.78*exp(17.27*(temperature - 273.15)/(temperature - 35.85));
    volScalarField pvap = relativeHumidity*psat;
    volScalarField pvap_surf = psat; // Assume saturation at droplet surface
    
    // Mass transfer coefficient
    volScalarField Dv = 2.5e-5; // Diffusivity of water vapor in air [m2/s]
    volScalarField km = Sh*Dv/diameter;
    
    // Evaporation rate [kg/m3/s]
    volScalarField Mw(0.018); // Molecular weight of water [kg/mol]
    volScalarField R(8.314); // Universal gas constant
    volScalarField evapRate = km*Mw*(pvap_surf - pvap)/(R*temperature);
    
    return evapRate;
}


Foam::volScalarField Foam::dropletPhysics::breakupTimeScale
(
    const volScalarField& diameter,
    const volScalarField& Urel,
    const volScalarField& rhoa,
    const volScalarField& rhop,
    const volScalarField& sigma
)
{
    // Calculate Weber number
    volScalarField We = WeberNumber(diameter, Urel, rhoa, sigma);
    
    // Breakup time scale (Pilch & Erdman 1987)
    volScalarField tBreakup = diameter*sqrt(rhop/rhoa)/mag(Urel);
    
    forAll(We, celli)
    {
        scalar We_val = We[celli];
        scalar factor = 1.0;
        
        if (We_val < 12)
        {
            // No breakup
            factor = GREAT;
        }
        else if (We_val < 18)
        {
            // Vibrational breakup
            factor = 6.0;
        }
        else if (We_val < 45)
        {
            // Bag breakup
            factor = 2.45*sqrt(We_val - 12);
        }
        else if (We_val < 100)
        {
            // Bag-and-stamen breakup
            factor = 14.1*pow(We_val - 12, -0.25);
        }
        else
        {
            // Shear breakup
            factor = 0.766*pow(We_val - 12, 0.25);
        }
        
        tBreakup[celli] *= factor;
    }
    
    return tBreakup;
}


Foam::volScalarField Foam::dropletPhysics::WeberNumber
(
    const volScalarField& diameter,
    const volScalarField& Urel,
    const volScalarField& rhoa,
    const volScalarField& sigma
)
{
    return rhoa*sqr(Urel)*diameter/sigma;
}


Foam::scalar Foam::dropletPhysics::criticalWeber(const scalar Oh)
{
    // Critical Weber number as function of Ohnesorge number
    // Based on Guildenbecher et al. (2009)
    if (Oh < 0.1)
    {
        return 12.0*(1.0 + 1.077*pow(Oh, 1.6));
    }
    else
    {
        return 12.0*(1.0 + 3.35*Oh);
    }
}


Foam::volScalarField Foam::dropletPhysics::OhnesorgeNumber
(
    const volScalarField& diameter,
    const volScalarField& mup,
    const volScalarField& rhop,
    const volScalarField& sigma
)
{
    return mup/sqrt(rhop*sigma*diameter);
}


Foam::scalar Foam::dropletPhysics::terminalVelocity(const scalar diameter)
{
    // Terminal velocity based on Gunn & Kinzer (1949) data
    // Simplified fit for diameters 0.1-5 mm
    scalar d_mm = diameter*1000.0;
    
    if (d_mm < 0.1)
    {
        return 0.27*d_mm;
    }
    else if (d_mm < 1.0)
    {
        return 4.03*pow(d_mm, 0.67);
    }
    else if (d_mm < 5.0)
    {
        return 9.17 - 0.99*d_mm + 0.106*d_mm*d_mm;
    }
    else
    {
        return 9.58; // Maximum for large drops
    }
}


Foam::volScalarField Foam::dropletPhysics::deformationDragFactor
(
    const volScalarField& We
)
{
    // Drag enhancement due to deformation
    // Based on Beard (1976)
    volScalarField dragFactor = 1.0 + 0.045*We;
    
    forAll(We, celli)
    {
        if (We[celli] > 8.0)
        {
            // Significant deformation
            dragFactor[celli] = 1.0 + 0.15*pow(We[celli], 0.687);
        }
    }
    
    return dragFactor;
}


Foam::scalar Foam::dropletPhysics::collisionEfficiency
(
    const scalar d1,
    const scalar d2,
    const scalar U1,
    const scalar U2
)
{
    // Collision efficiency for two drops
    // Based on Beard & Ochs (1984)
    scalar p = d2/d1; // Size ratio (smaller/larger)
    scalar St = 2.0*p*p*mag(U1 - U2)*d2/(9.0*1.8e-5*d1); // Stokes number
    
    if (St < 0.083)
    {
        return 0.0;
    }
    else if (St < 1.0)
    {
        return 4.5*St*St;
    }
    else
    {
        scalar StCrit = 1.0 + 0.75*log(2.0*p);
        if (St > StCrit)
        {
            return pow(St - StCrit + 1.0, -2.0);
        }
        else
        {
            return 1.0;
        }
    }
}


Foam::scalar Foam::dropletPhysics::coalescenceEfficiency
(
    const scalar d1,
    const scalar d2,
    const scalar Weber
)
{
    // Coalescence efficiency
    // Based on Low & List (1982)
    scalar p = min(d1, d2)/max(d1, d2);
    
    if (Weber < 1.0)
    {
        return 1.0; // Always coalesce at low Weber
    }
    else if (Weber < 7.0)
    {
        return exp(-0.5*pow((Weber - 1.0)/3.0, 2));
    }
    else
    {
        return 0.0; // Bounce or breakup
    }
}


Foam::volScalarField Foam::dropletPhysics::evaporationHeatSource
(
    const volScalarField& evapRate,
    const volScalarField& temperature
)
{
    // Latent heat of vaporization [J/kg]
    scalar Lv = 2.5e6 - 2370.0*(temperature.average().value() - 273.15);
    
    // Heat source due to evaporation [W/m3]
    return -evapRate*Lv;
}


// ************************************************************************* //