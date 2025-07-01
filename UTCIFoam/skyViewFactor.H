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
    Sky View Factor (SVF) calculation using ray tracing method
    
    SVF represents the fraction of sky visible from a point, ranging from 0 to 1:
    - 1.0 = completely open (unobstructed view of sky)
    - 0.0 = completely obstructed (no sky visible)
    
    Method: Ray casting from each cell center in multiple directions to detect
    sky dome visibility

\*---------------------------------------------------------------------------*/

#ifndef skyViewFactor_H
#define skyViewFactor_H

#include "volFields.H"
#include "meshSearch.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate SVF using ray tracing
volScalarField calculateSkyViewFactor
(
    const fvMesh& mesh,
    const label nRays = 100,        // Number of rays for sampling
    const scalar maxDistance = 1000 // Maximum ray distance [m]
)
{
    Info<< "Calculating Sky View Factor using " << nRays << " rays per cell..." << endl;
    
    volScalarField SVF
    (
        IOobject
        (
            "SVF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SVF", dimless, 0.0)
    );
    
    // Create mesh search engine
    meshSearch searchEngine(mesh);
    
    // Random number generator for ray directions
    Random rndGen(1234567);
    
    // Get mesh bounds to determine domain top
    const boundBox& meshBb = mesh.bounds();
    const scalar domainTop = meshBb.max().z() + maxDistance;
    
    // Process each cell
    forAll(SVF, cellI)
    {
        const point& cellCenter = mesh.C()[cellI];
        label nHits = 0;
        
        // Generate rays in hemisphere above the point
        for (label rayI = 0; rayI < nRays; rayI++)
        {
            // Generate random direction in upper hemisphere
            // Using stratified random sampling for better distribution
            
            // Azimuth angle (0 to 2*pi)
            scalar phi = 2.0 * constant::mathematical::pi * rndGen.sample01<scalar>();
            
            // Elevation angle (0 to pi/2) with cosine weighting
            // This ensures uniform distribution on hemisphere
            scalar cosTheta = rndGen.sample01<scalar>();  // cos(theta) from 0 to 1
            scalar sinTheta = Foam::sqrt(1.0 - cosTheta*cosTheta);
            
            // Convert to Cartesian direction
            vector rayDir
            (
                sinTheta * Foam::cos(phi),
                sinTheta * Foam::sin(phi),
                cosTheta  // Pointing upward
            );
            
            // Ray end point
            point rayEnd = cellCenter + maxDistance * rayDir;
            
            // Check if ray hits any face before reaching domain top
            pointIndexHit hitInfo = searchEngine.intersection(cellCenter, rayEnd);
            
            if (!hitInfo.hit())
            {
                // Ray reached sky without obstruction
                nHits++;
            }
            else
            {
                // Check if hit point is at domain boundary (sky)
                const point& hitPoint = hitInfo.hitPoint();
                if (hitPoint.z() >= domainTop * 0.95)  // Near top boundary
                {
                    nHits++;
                }
            }
        }
        
        // Calculate SVF as ratio of unobstructed rays
        SVF[cellI] = scalar(nHits) / scalar(nRays);
    }
    
    // Apply boundary conditions
    SVF.correctBoundaryConditions();
    
    Info<< "Sky View Factor calculation completed." << endl;
    Info<< "SVF range: " << min(SVF).value() << " to " << max(SVF).value() << endl;
    
    return SVF;
}

// Simplified SVF calculation based on vertical obstruction angle
volScalarField calculateSimplifiedSVF
(
    const fvMesh& mesh,
    const scalar defaultSVF = 0.7  // Default for urban environments
)
{
    Info<< "Using simplified Sky View Factor calculation..." << endl;
    
    volScalarField SVF
    (
        IOobject
        (
            "SVF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SVF", dimless, defaultSVF)
    );
    
    // Get mesh bounds
    const boundBox& meshBb = mesh.bounds();
    const scalar groundLevel = meshBb.min().z();
    const scalar maxHeight = meshBb.max().z();
    const scalar domainHeight = maxHeight - groundLevel;
    
    // Simple height-based estimation
    forAll(SVF, cellI)
    {
        const scalar cellHeight = mesh.C()[cellI].z() - groundLevel;
        const scalar relativeHeight = cellHeight / domainHeight;
        
        // Higher cells have better sky view
        // Ground level: SVF ≈ 0.5 (urban canyon)
        // Top level: SVF ≈ 1.0 (unobstructed)
        SVF[cellI] = 0.5 + 0.5 * relativeHeight;
        
        // Clamp to valid range
        SVF[cellI] = max(0.0, min(1.0, SVF[cellI]));
    }
    
    return SVF;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //