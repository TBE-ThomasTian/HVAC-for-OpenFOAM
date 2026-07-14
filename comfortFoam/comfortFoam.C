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

Application
    comfortFoam

Description
    This tool calculates thermal comfort parameters according to
    ISO 7730:2025 and DIN EN ISO 7730:2025-12:
    - PMV (Predicted Mean Vote), model applicability: -2 to +2
    - PPD (Predicted Percentage of Dissatisfied), Range: 0 to 100%
    - DR  (Draught Rating), Range: 0 to 100%
    - TOp (Operative Temperature)
    - Local discomfort due to vertical temperature difference, floor
      temperature and radiant temperature asymmetry when inputs are available
    - ISO 7730 whole-body and overall categories
    - Mean radiant temperature
    - Turbulent intensity %

    Analysis can be limited to specific regions:
    - Use -cellSet <name> (or -setFields <name>) for a specific cellSet
    - Use -cellZone <name> for a specific cellZone
    - Without options: analyzes entire mesh

    For humidity calculations, run: buoyantHumiditySimpleFoam / buoyantHumidityPimpleFoam
    
    Supported humidity fields:
    - thermo:relHum (from buoyantHumiditySimpleFoam, range [0,1])
    - thermoRelHum (legacy naming, range [0,1])
    - relHum (legacy field, range [0,100] or [0,1])
    - Default constant value if no field available
    
    Supported radiation fields (in order of preference):
    - MRT: Mean radiant temperature (K)
    - G: Incident radiation field (W/m^2) from radiation models
    - qr: Radiative heat flux field (W/m^2) from some radiation models
    - IDefault: Default radiation intensity (W/m^2/sr) from DOM models
    - Fallback: Area-weighted wall temperature calculation

    Usage examples:
    comfortFoam                           # Analyze entire mesh
    comfortFoam -cellSet occupancyZone    # Analyze only specified cellSet
    comfortFoam -setFields occupancyZone  # Alias for -cellSet
    comfortFoam -cellZone roomA           # Analyze only specified cellZone

Background
    ISO 7730:2025; DIN EN ISO 7730:2025-12

Authors
    Thomas Tian
    Tobias Holzmann
    Manuel Scheu

Version
    5.0

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "wallFvPatch.H"
#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "cellSet.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "atmBoundaryLayerInletVelocityFvPatchVectorField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Constants * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace ComfortConstants
    {
        // Physical constants
        const scalar stefanBoltzmannConstant = 5.67e-8;  // W/(m^2*K^4)
        const scalar baseMetabolicRate = 58.15;          // W/m^2
        const scalar physicalKelvinOffset = 273.15;      // K at 0 degrees C
        const scalar isoTemperatureOffset = 273.0;       // Annex D program
        
        // Heat loss coefficients (from ISO 7730)
        const scalar skinDiffusionCoeff = 3.05e-3;
        const scalar thermalCoeff = 6.99;
        const scalar basePressure = 5733;                // Pa
        const scalar respirationCoeff1 = 1.7e-5;
        const scalar respirationCoeff2 = 0.0014;
        const scalar respirationTemp = 34.0;             // degrees C
        const scalar respirationHumidity = 5867;         // Pa
        const scalar radiationCoeff = 3.96;
        const scalar thermalSensCoeff1 = 0.303;
        const scalar thermalSensCoeff2 = 0.036;  // ISO 7730 exponential coefficient
        const scalar thermalSensCoeff3 = 0.028;  // ISO 7730 additive constant
        
        // Clothing coefficients
        const scalar clothingFactor1 = 1.29;
        const scalar clothingFactor2 = 1.05;
        const scalar clothingFactor3 = 0.645;
        const scalar clothingThreshold = 0.078;
        
        // Convection coefficients
        const scalar forcedConvectionCoeff = 12.1;
        const scalar naturalConvectionCoeff = 2.38;
        const scalar naturalConvectionExp = 0.25;
        
        // Draft calculation coefficients
        const scalar draftBaseTemp = 34.0;               // degrees C
        const scalar draftVelThreshold = 0.05;           // m/s
        const scalar draftVelExp = 0.62;
        const scalar draftTurbCoeff = 0.37;
        const scalar draftConstant = 3.14;
        
        // Iteration parameters
        const label maxClothingIterations = 150;
        const scalar clothingConvergenceTol = 0.00015;  // ISO 7730 epsilon from BASIC line 340
        
        // PPD calculation coefficients
        const scalar ppdCoeff1 = 0.03353;
        const scalar ppdCoeff2 = 0.2179;
        const scalar ppdBase = 95.0;
        
        // Operative temperature correction factors
        const scalar opTempFactor1 = 0.5;  // v < 0.2 m/s
        const scalar opTempFactor2 = 0.6;  // 0.2 <= v <= 0.6 m/s
        const scalar opTempFactor3 = 0.7;  // v > 0.6 m/s
        const scalar velThreshold1 = 0.2;  // m/s
        const scalar velThreshold2 = 0.6;  // m/s

        // ISO 7730 PMV applicability limits
        const scalar minAirTemperature = 10.0;       // degrees C
        const scalar maxAirTemperature = 30.0;       // degrees C
        const scalar minRadiantTemperature = 10.0;   // degrees C
        const scalar maxRadiantTemperature = 40.0;   // degrees C
        const scalar minRelativeAirVelocity = 0.0;   // m/s
        const scalar maxRelativeAirVelocity = 1.0;   // m/s
        const scalar minWaterVapourPressure = 0.0;   // Pa
        const scalar maxWaterVapourPressure = 2700.0;// Pa
        const scalar minApplicablePMV = -2.0;
        const scalar maxApplicablePMV = 2.0;

        // ISO 7730 draught-model applicability limits
        const scalar minDraftTemperature = 20.0;     // degrees C
        const scalar maxDraftTemperature = 26.0;     // degrees C
        const scalar maxDraftVelocity = 0.5;         // m/s, exclusive
        const scalar minDraftTurbulence = 10.0;      // percent
        const scalar maxDraftTurbulence = 60.0;      // percent

        // ISO 9920 / ISO 7730 Annex C dynamic-clothing constants
        const scalar boundaryAirLayerInsulation = 0.111; // m2*K/W
    }

// * * * * * * * * * * * * * * * Enums * * * * * * * * * * * * * * * * * * //

enum class ComfortCategory
{
    Incomplete = -1,
    None = 0,
    I = 1,
    II = 2,
    III = 3,
    IV = 4
};

enum class AnalysisRegionType
{
    EntireMesh,
    CellSet,
    CellZone
};

// Convert comfort category to an English description.
const char* comfortCategoryToString(ComfortCategory category)
{
    switch (category)
    {
        case ComfortCategory::I: return "Category I";
        case ComfortCategory::II: return "Category II";
        case ComfortCategory::III: return "Category III";
        case ComfortCategory::IV: return "Category IV";
        case ComfortCategory::Incomplete:
            return "Not determined (local-discomfort input missing)";
        default: return "No category";
    }
}

Foam::scalar comfortCategoryValue(ComfortCategory category)
{
    return static_cast<label>(category);
}

// Calculate volume-weighted averages for specific cells
void calculateVolumeWeightedAverages
(
    const fvMesh& mesh,
    const labelList& cellsToAnalyze,
    scalar& avgVelocity,
    scalar& avgTemperature,
    scalar& totalVolume
)
{
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    
    scalar weightedVelocityMagnitude(0);
    scalar weightedTemperature(0);
    totalVolume = 0;
    
    forAll(cellsToAnalyze, i)
    {
        label cellI = cellsToAnalyze[i];
        scalar cellVolume = mesh.V()[cellI];
        
        weightedVelocityMagnitude += mag(U[cellI]) * cellVolume;
        weightedTemperature += T[cellI] * cellVolume;
        totalVolume += cellVolume;
    }
    
    // Reduce for parallel operation
    reduce(weightedVelocityMagnitude, sumOp<scalar>());
    reduce(weightedTemperature, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    
    avgVelocity = weightedVelocityMagnitude / totalVolume;
    avgTemperature = weightedTemperature / totalVolume;
}
}

// * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

// Validate input parameters
void validateInputParameters(scalar met, scalar clo, scalar wme, scalar rh)
{
    if (met < 0.8 || met > 4.0)
    {
        FatalErrorInFunction
            << "Metabolic rate out of valid range (0.8-4.0 met): " << met
            << abort(FatalError);
    }
    
    if (clo < 0 || clo > 2.0)
    {
        FatalErrorInFunction 
            << "Clothing insulation out of valid range (0-2.0 clo): " << clo
            << abort(FatalError);
    }
    
    if (wme < 0 || wme > met)
    {
        FatalErrorInFunction
            << "External work out of valid range (0-" << met << " met): " << wme
            << abort(FatalError);
    }
    
    if (rh < 0 || rh > 100)
    {
        FatalErrorInFunction
            << "Relative humidity out of valid range (0-100%): " << rh
            << abort(FatalError);
    }
}

// Calculate volume-weighted average velocity
Foam::vector calculateAverageVelocity(const fvMesh& mesh)
{
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    vector avgU(Zero);
    scalar totalVolume(0);

    forAll(mesh.cells(), cellI)
    {
        avgU += U[cellI] * mesh.V()[cellI];
        totalVolume += mesh.V()[cellI];
    }

    return avgU / totalVolume;
}

// Calculate volume-weighted average temperature
Foam::scalar calculateAverageTemperature(const fvMesh& mesh)
{
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    scalar avgT(0);
    scalar totalVolume(0);

    forAll(mesh.cells(), cellI)
    {
        avgT += T[cellI] * mesh.V()[cellI];
        totalVolume += mesh.V()[cellI];
    }

    return avgT / totalVolume;
}

// Approximate mean radiant temperature from wall temperatures. This fallback
// uses fourth-power area weighting. A measured or view-factor-based MRT field
// is preferred because this approximation is not occupant-directional.
Foam::scalar calculateRadiationTemperature
(
    const fvMesh& mesh,
    const fvPatchList& patches
)
{
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    scalar weightedTemperature4(0);
    scalar totalArea(0);

    forAll(patches, patchI)
    {
        const label curPatch = patches[patchI].index();

        if (isType<wallFvPatch>(patches[patchI]))
        {
            scalar patchArea = gSum(mesh.magSf().boundaryField()[curPatch]);

            if (patchArea > SMALL)
            {
                weightedTemperature4 += gSum
                (
                    mesh.magSf().boundaryField()[curPatch]
                  * Foam::pow(T.boundaryField()[curPatch], 4)
                );

                totalArea += patchArea;
            }
        }
    }

    if (totalArea <= SMALL)
    {
        WarningInFunction
            << "No wall area is available for the MRT fallback. "
            << "Using the volume-average air temperature." << endl;
        return
            gAverage(T.internalField())
          - ComfortConstants::physicalKelvinOffset;
    }

    return
        Foam::pow(weightedTemperature4/totalArea, 0.25)
      - ComfortConstants::physicalKelvinOffset;
}

// Calculate saturation water-vapour pressure term used by Annex D.
Foam::scalar calculateSaturationPressureTerm(scalar airTemperatureC)
{
    return Foam::exp(16.6536 - (4030.183/(airTemperatureC + 235.0)));
}

// Calculate water-vapour partial pressure [Pa] from RH [%].
Foam::scalar calculateWaterVapourPressure
(
    scalar airTemperatureC,
    scalar relativeHumidity
)
{
    return
        relativeHumidity
       *10.0
       *calculateSaturationPressureTerm(airTemperatureC);
}

// Calculate RH [%] for reporting when water-vapour pressure is supplied.
Foam::scalar calculateRelativeHumidity
(
    scalar airTemperatureC,
    scalar vapourPressure
)
{
    const scalar denominator =
        10.0*calculateSaturationPressureTerm(airTemperatureC);
    return denominator > SMALL ? 100.0*vapourPressure/denominator : 0.0;
}

// Calculate turbulence intensity as a fraction (0.4 means 40 percent).
Foam::scalar calculateTurbulentIntensityFraction
(
    scalar velocity,
    scalar k,
    bool turbulenceAvailable
)
{
    if (turbulenceAvailable && velocity > SMALL)
    {
        return Foam::sqrt((2.0/3.0)*Foam::max(k, 0.0))/velocity;
    }
    return 0.4;
}

// Calculate draught rating (DR)
Foam::scalar calculateDraughtRating
(
    scalar airTemperatureC,
    scalar velocity,
    scalar turbulentIntensityFraction
)
{
    using namespace ComfortConstants;
    
    const scalar modelVelocity = Foam::max(velocity, draftVelThreshold);
    const scalar turbulentIntensityPercent =
        100.0*turbulentIntensityFraction;

    const scalar dr = (draftBaseTemp - airTemperatureC)
                    * Foam::pow(modelVelocity - draftVelThreshold, draftVelExp)
                    *
                      (
                          draftTurbCoeff*modelVelocity
                         *turbulentIntensityPercent
                        + draftConstant
                      );

    return Foam::max(0.0, Foam::min(100.0, dr));
}

bool isDraughtModelApplicable
(
    scalar airTemperatureC,
    scalar velocity,
    scalar turbulentIntensityFraction,
    scalar met
)
{
    using namespace ComfortConstants;

    const scalar turbulencePercent = 100.0*turbulentIntensityFraction;

    return
        airTemperatureC >= minDraftTemperature
     && airTemperatureC <= maxDraftTemperature
     && velocity < maxDraftVelocity
     && turbulencePercent >= minDraftTurbulence
     && turbulencePercent <= maxDraftTurbulence
     && met <= 1.2;
}

Foam::scalar calculateClothingAreaFactor(scalar insulationSI)
{
    using namespace ComfortConstants;

    return insulationSI <= clothingThreshold
        ? 1.0 + clothingFactor1*insulationSI
        : clothingFactor2 + clothingFactor3*insulationSI;
}

Foam::scalar calculateInitialClothingSurfaceTemperature
(
    scalar airTemp,
    scalar airTempC,
    scalar insulationSI
)
{
    // ISO 7730:2025 Annex D, BASIC line 250.
    return airTemp
         + (35.5 - airTempC)/(3.5*(6.45*insulationSI + 0.1));
}

Foam::scalar calculateWalkingSpeed
(
    scalar metabolicRate,
    scalar configuredWalkingSpeed
)
{
    if (configuredWalkingSpeed >= 0.0)
    {
        return Foam::min(configuredWalkingSpeed, 1.2);
    }

    // ISO 7730:2025 Annex C / ISO 9920 activity-based estimate.
    return Foam::min
    (
        0.7,
        Foam::max(0.0, 0.0052*(metabolicRate - 58.0))
    );
}

Foam::scalar calculateResultantClothingInsulation
(
    scalar staticClothingInsulation,
    scalar relativeAirVelocity,
    scalar walkingSpeed,
    bool& correctionApplicable
)
{
    using namespace ComfortConstants;

    correctionApplicable =
        staticClothingInsulation >= 0.0
     && staticClothingInsulation < 1.4;

    if (!correctionApplicable)
    {
        return staticClothingInsulation;
    }

    const scalar limitedVelocity =
        Foam::min(Foam::max(relativeAirVelocity, 0.15), 3.5);
    const scalar limitedWalkingSpeed =
        Foam::min(Foam::max(walkingSpeed, 0.0), 1.2);
    const scalar velocityOffset = limitedVelocity - 0.15;

    // ISO 9920:2007, equations 32 and 33. The 0.044 coefficient is
    // intentional. The 0.44 value in ISO 7730 Annex C is an apparent
    // typographical discrepancy that suppresses the expected wind correction.
    const scalar correctionTotal = Foam::min
    (
        1.0,
        Foam::exp
        (
            -0.281*velocityOffset
          + 0.044*Foam::sqr(velocityOffset)
          - 0.492*limitedWalkingSpeed
          + 0.176*Foam::sqr(limitedWalkingSpeed)
        )
    );
    const scalar correctionAirLayer = Foam::min
    (
        1.0,
        Foam::exp
        (
            -0.533*velocityOffset
          + 0.069*Foam::sqr(velocityOffset)
          - 0.462*limitedWalkingSpeed
          + 0.201*Foam::sqr(limitedWalkingSpeed)
        )
    );

    const scalar targetInsulationSI = 0.155*staticClothingInsulation;
    const scalar targetAreaFactor =
        calculateClothingAreaFactor(targetInsulationSI);
    const scalar resultantAirLayer =
        boundaryAirLayerInsulation*correctionAirLayer;

    scalar resultantTotalInsulation = 0.0;

    if (staticClothingInsulation <= 0.6)
    {
        // ISO 9920 equation 34 interpolates between the nude and 0.6-clo
        // endpoint totals, not between corrections evaluated at target clo.
        const scalar endpointClo = 0.6;
        const scalar endpointInsulationSI = 0.155*endpointClo;
        const scalar endpointAreaFactor =
            calculateClothingAreaFactor(endpointInsulationSI);
        const scalar endpointTotalInsulation =
            endpointInsulationSI
          + boundaryAirLayerInsulation/endpointAreaFactor;
        const scalar resultantNudeTotal = resultantAirLayer;
        const scalar resultantDressedTotal =
            correctionTotal*endpointTotalInsulation;

        resultantTotalInsulation =
        (
            (endpointClo - staticClothingInsulation)*resultantNudeTotal
          + staticClothingInsulation*resultantDressedTotal
        )/endpointClo;
    }
    else
    {
        const scalar totalInsulation =
            targetInsulationSI
          + boundaryAirLayerInsulation/targetAreaFactor;
        resultantTotalInsulation = correctionTotal*totalInsulation;
    }

    const scalar resultantClothingSI = Foam::max
    (
        0.0,
        resultantTotalInsulation - resultantAirLayer/targetAreaFactor
    );

    return resultantClothingSI/0.155;
}

// Solve for clothing surface temperature iteratively
// Returns Tuple2<scalar> with first = tcl (Celsius), second = xn (Kelvin/100)
Foam::Tuple2<scalar> calculateClothingSurfaceTemperature
(
    scalar airTempC,
    scalar velocity,
    scalar icl,
    scalar fcl,
    scalar radiationTemp,
    scalar metabolicRate,
    scalar externalWork,
    bool& converged
)
{
    using namespace ComfortConstants;

    // Annex D intentionally uses 273.0 internally. OpenFOAM field values are
    // converted from physical kelvin to degrees Celsius before this function.
    const scalar isoAirTemperature = airTempC + isoTemperatureOffset;
    const scalar tcla = calculateInitialClothingSurfaceTemperature
    (
        isoAirTemperature,
        airTempC,
        icl
    );
    // ISO 7730: iteration uses Kelvin divided by 100
    scalar xn = tcla / 100.0;
    scalar xf = xn;
    
    scalar p1 = icl * fcl;
    scalar p2 = p1 * radiationCoeff;
    scalar p3 = p1 * 100.0;
    scalar p4 = p1 * isoAirTemperature;
    scalar p5 = 308.7 - 0.028 * (metabolicRate - externalWork)
              + p2
               *Foam::pow
                (
                    (radiationTemp + isoTemperatureOffset)/100.0,
                    4
                );
    
    label iterationCount = 0;
    converged = false;
    
    
    do
    {
        iterationCount++;
        xf = (xf + xn) / 2.0;  // ISO 7730 BASIC line 350
        
        scalar hcf = forcedConvectionCoeff * Foam::sqrt(velocity);
        scalar hcn = naturalConvectionCoeff
            *Foam::pow
             (
                 mag(100.0*xf - isoAirTemperature),
                 naturalConvectionExp
             );
        scalar hc = Foam::max(hcf, hcn);
        
        xn = (p5 + p4 * hc - p2 * Foam::pow(xf, 4.0)) / (100.0 + p3 * hc);
        
        
        if (iterationCount > maxClothingIterations)
        {
            break;
        }

        converged = mag(xn - xf) <= clothingConvergenceTol;

    } while (!converged);
    
    // Return both tcl (Celsius) and xn (Kelvin/100)
    return Tuple2<scalar>(100.0*xn - isoTemperatureOffset, xn);
}

// Calculate PMV (Predicted Mean Vote)
Foam::scalar calculatePMV
(
    scalar airTempC,
    scalar velocity,
    scalar waterVapourPressure,
    scalar radiationTemp,
    scalar metabolicRate,
    scalar clothingInsulation,
    scalar externalWork,
    bool& converged
)
{
    using namespace ComfortConstants;

    scalar icl = 0.155 * clothingInsulation;
    scalar fcl = calculateClothingAreaFactor(icl);
    
    Tuple2<scalar> tclResult = calculateClothingSurfaceTemperature
    (
        airTempC,
        velocity,
        icl,
        fcl,
        radiationTemp,
        metabolicRate,
        externalWork,
        converged
    );
    scalar tcl = tclResult.first();  // Temperature in Celsius
    scalar xn = tclResult.second();   // Kelvin/100
    
    // Calculate heat losses
    scalar hl1 = skinDiffusionCoeff
        *
        (
            basePressure
          - thermalCoeff*(metabolicRate - externalWork)
          - waterVapourPressure
        );
    
    scalar hl2 = 0.0;
    if ((metabolicRate - externalWork) > baseMetabolicRate)
    {
        hl2 = 0.42 * ((metabolicRate - externalWork) - baseMetabolicRate);
    }
    
    scalar hl3 = respirationCoeff1*metabolicRate
        *(respirationHumidity - waterVapourPressure);
    scalar hl4 = respirationCoeff2 * metabolicRate * (respirationTemp - airTempC);
    
    // Use xn directly for HL5 calculation (like ISO 7730 BASIC line 480)
    scalar hl5 = radiationCoeff * fcl * 
                (
                    Foam::pow(xn, 4)
                  - Foam::pow
                    (
                        (radiationTemp + isoTemperatureOffset)/100.0,
                        4
                    )
                );
    
    scalar hcf = forcedConvectionCoeff * Foam::sqrt(velocity);
    scalar hcn = naturalConvectionCoeff * Foam::pow(mag(tcl - airTempC), naturalConvectionExp);
    scalar hc = Foam::max(hcf, hcn);
    
    scalar hl6 = fcl * hc * (tcl - airTempC);
    
    // ISO 7730 formula: TS = 0.303 * exp(-0.036 * M) + 0.028
    scalar ts = thermalSensCoeff1 * Foam::exp(-thermalSensCoeff2 * metabolicRate) + thermalSensCoeff3;
    
    return ts * ((metabolicRate - externalWork) - hl1 - hl2 - hl3 - hl4 - hl5 - hl6);
}

// Calculate PPD (Predicted Percentage of Dissatisfied)
Foam::scalar calculatePPD(scalar pmv)
{
    using namespace ComfortConstants;
    
    return 100.0 - ppdBase * Foam::exp(-ppdCoeff1 * Foam::pow(pmv, 4) - ppdCoeff2 * Foam::pow(pmv, 2));
}

bool isPMVModelApplicable
(
    scalar airTemperatureC,
    scalar radiantTemperature,
    scalar relativeAirVelocity,
    scalar waterVapourPressure,
    scalar clothingInsulation,
    scalar pmv,
    bool clothingIterationConverged,
    bool clothingInputCompliant
)
{
    using namespace ComfortConstants;

    return
        clothingIterationConverged
     && clothingInputCompliant
     && airTemperatureC >= minAirTemperature
     && airTemperatureC <= maxAirTemperature
     && radiantTemperature >= minRadiantTemperature
     && radiantTemperature <= maxRadiantTemperature
     && relativeAirVelocity >= minRelativeAirVelocity
     && relativeAirVelocity <= maxRelativeAirVelocity
     && waterVapourPressure >= minWaterVapourPressure
     && waterVapourPressure <= maxWaterVapourPressure
     && clothingInsulation >= 0.0
     && clothingInsulation <= 2.0
     && pmv >= minApplicablePMV
     && pmv <= maxApplicablePMV;
}

Foam::scalar calculateVerticalTemperatureDissatisfaction
(
    scalar verticalTemperatureDifference
)
{
    return 100.0
        /
        (
            1.0
          + Foam::exp(5.76 - 0.856*verticalTemperatureDifference)
        );
}

Foam::scalar calculateFloorTemperatureDissatisfaction
(
    scalar floorTemperatureC
)
{
    const scalar pd = 100.0
        - 94.0
         *Foam::exp
          (
              -1.387
            + 0.118*floorTemperatureC
            - 0.0025*Foam::sqr(floorTemperatureC)
          );
    return Foam::max(0.0, Foam::min(100.0, pd));
}

bool radiantAsymmetryConstants
(
    const word& asymmetryType,
    scalar& k1,
    scalar& k2,
    scalar& k3,
    scalar& upperLimit
)
{
    if (asymmetryType == "warmCeiling")
    {
        k1 = 2.94;
        k2 = 0.166;
        k3 = 5.5;
        upperLimit = 23.0;
    }
    else if (asymmetryType == "coolWall")
    {
        k1 = 5.89;
        k2 = 0.297;
        k3 = 1.0;
        upperLimit = 15.0;
    }
    else if (asymmetryType == "coolCeiling")
    {
        k1 = 5.19;
        k2 = 0.173;
        k3 = 1.0;
        upperLimit = 15.0;
    }
    else if (asymmetryType == "warmWall")
    {
        k1 = 3.41;
        k2 = 0.044;
        k3 = 3.5;
        upperLimit = 35.0;
    }
    else
    {
        return false;
    }

    return true;
}

Foam::scalar calculateRadiantAsymmetryDissatisfaction
(
    scalar radiantTemperatureAsymmetry,
    scalar k1,
    scalar k2,
    scalar k3
)
{
    return Foam::max
    (
        0.0,
        Foam::min
        (
            100.0,
            100.0
           /
            (
                1.0
              + Foam::exp(k1 - k2*radiantTemperatureAsymmetry)
            )
          - k3
        )
    );
}

// Calculate operative temperature
Foam::scalar calculateOperativeTemperature
(
    scalar airTemperatureC,
    scalar radiationTemperatureC,
    scalar velocity
)
{
    using namespace ComfortConstants;
    
    scalar correctionFactor;
    
    if (velocity < velThreshold1)
    {
        correctionFactor = opTempFactor1;
    }
    else if (velocity <= velThreshold2)  // Fixed logical error: AND instead of OR
    {
        correctionFactor = opTempFactor2;
    }
    else
    {
        correctionFactor = opTempFactor3;
    }
    
    return
        correctionFactor*airTemperatureC
      + (1.0 - correctionFactor)*radiationTemperatureC
      + physicalKelvinOffset;
}

// Annex A categories are informative. Strict inequalities follow Table A.1.
ComfortCategory analyzeWholeBodyCategory(scalar pmv, scalar ppd)
{
    if (pmv > -0.2 && pmv < 0.2 && ppd < 6.0)
    {
        return ComfortCategory::I;
    }
    if (pmv > -0.5 && pmv < 0.5 && ppd < 10.0)
    {
        return ComfortCategory::II;
    }
    if (pmv > -0.7 && pmv < 0.7 && ppd < 15.0)
    {
        return ComfortCategory::III;
    }
    if (pmv > -1.0 && pmv < 1.0 && ppd < 25.0)
    {
        return ComfortCategory::IV;
    }

    return ComfortCategory::None;
}

ComfortCategory analyzeOverallComfortCategory
(
    scalar pmv,
    scalar ppd,
    scalar dr,
    scalar verticalPD,
    scalar floorPD,
    scalar radiantAsymmetryPD,
    bool localInputsAvailable
)
{
    const ComfortCategory wholeBody = analyzeWholeBodyCategory(pmv, ppd);

    if (!localInputsAvailable)
    {
        return wholeBody == ComfortCategory::IV
            ? ComfortCategory::IV
            : ComfortCategory::Incomplete;
    }

    if
    (
        pmv > -0.2 && pmv < 0.2
     && ppd < 6.0
     && dr < 10.0
     && verticalPD < 3.0
     && floorPD < 10.0
     && radiantAsymmetryPD < 5.0
    )
    {
        return ComfortCategory::I;
    }

    if
    (
        pmv > -0.5 && pmv < 0.5
     && ppd < 10.0
     && dr < 20.0
     && verticalPD < 5.0
     && floorPD < 10.0
     && radiantAsymmetryPD < 5.0
    )
    {
        return ComfortCategory::II;
    }

    if
    (
        pmv > -0.7 && pmv < 0.7
     && ppd < 15.0
     && dr < 30.0
     && verticalPD < 10.0
     && floorPD < 15.0
     && radiantAsymmetryPD < 10.0
    )
    {
        return ComfortCategory::III;
    }

    return wholeBody == ComfortCategory::None
        ? ComfortCategory::None
        : ComfortCategory::IV;
}

// Normalize selection labels to valid local cell labels.
// In parallel, this also handles the case where a set/zone stores global labels.
Foam::labelList normalizeSelectedCells
(
    const fvMesh& mesh,
    const labelUList& rawCells,
    const word& selectionKind,
    const word& selectionName
)
{
    const label nLocalCells = mesh.nCells();
    bool hasOutOfLocalRange = false;

    forAll(rawCells, i)
    {
        const label cellI = rawCells[i];
        if (cellI < 0 || cellI >= nLocalCells)
        {
            hasOutOfLocalRange = true;
            break;
        }
    }

    if (Pstream::parRun())
    {
        reduce(hasOutOfLocalRange, orOp<bool>());
    }

    DynamicList<label> validCells(rawCells.size());
    label nDiscarded = 0;

    if (Pstream::parRun() && hasOutOfLocalRange)
    {
        const globalIndex& globalCellAddr = mesh.globalData().globalMeshCellAddr();

        forAll(rawCells, i)
        {
            const label globalCellI = rawCells[i];
            if (globalCellAddr.isLocal(globalCellI))
            {
                validCells.append(globalCellAddr.toLocal(globalCellI));
            }
            else
            {
                ++nDiscarded;
            }
        }

        const label totalDiscarded = returnReduce(nDiscarded, sumOp<label>());
        if (totalDiscarded > 0)
        {
            Info<< "Converted global labels for " << selectionKind << " "
                << selectionName << " and ignored " << totalDiscarded
                << " non-local entries across all processors." << endl;
        }
    }
    else
    {
        forAll(rawCells, i)
        {
            const label localCellI = rawCells[i];
            if (localCellI >= 0 && localCellI < nLocalCells)
            {
                validCells.append(localCellI);
            }
            else
            {
                ++nDiscarded;
            }
        }

        const label totalDiscarded = returnReduce(nDiscarded, sumOp<label>());
        if (totalDiscarded > 0)
        {
            WarningInFunction
                << "Ignoring " << totalDiscarded << " invalid labels in "
                << selectionKind << " " << selectionName << endl;
        }
    }

    labelList cellsToAnalyze(validCells.size());
    forAll(validCells, i)
    {
        cellsToAnalyze[i] = validCells[i];
    }

    // Defensive duplicate removal to avoid double-counting volume.
    if (!cellsToAnalyze.empty())
    {
        sort(cellsToAnalyze);
        label writeI = 0;
        for (label readI = 1; readI < cellsToAnalyze.size(); ++readI)
        {
            if (cellsToAnalyze[readI] != cellsToAnalyze[writeI])
            {
                cellsToAnalyze[++writeI] = cellsToAnalyze[readI];
            }
        }
        cellsToAnalyze.setSize(writeI + 1);
    }

    return cellsToAnalyze;
}

// Get list of cells to analyze (from CLI/dictionary cellSet/cellZone or entire mesh)
Foam::labelList getCellsToAnalyze
(
    const fvMesh& mesh,
    const argList& args,
    const dictionary& comfortFoamDict,
    AnalysisRegionType& selectionType,
    word& selectionName
)
{
    labelList cellsToAnalyze;
    selectionType = AnalysisRegionType::EntireMesh;
    selectionName.clear();

    const bool cliCellSetFound = args.found("cellSet");
    const bool cliSetFieldsFound = args.found("setFields");
    const bool cliCellZoneFound = args.found("cellZone");

    word cliCellSetName;
    if (cliCellSetFound)
    {
        cliCellSetName = args.get<word>("cellSet");
    }
    if (cliSetFieldsFound)
    {
        const word setFieldsName = args.get<word>("setFields");
        if (cliCellSetFound && setFieldsName != cliCellSetName)
        {
            FatalErrorInFunction
                << "Both -cellSet and -setFields were specified with different names: "
                << cliCellSetName << " vs " << setFieldsName << nl
                << "Use only one option (or use the same name)."
                << abort(FatalError);
        }
        cliCellSetName = setFieldsName;
    }

    const bool dictCellSetFound = comfortFoamDict.found("cellSet");
    const bool dictSetFieldsFound = comfortFoamDict.found("setFields");
    const bool dictCellZoneFound = comfortFoamDict.found("cellZone");

    word dictCellSetName;
    if (dictCellSetFound)
    {
        dictCellSetName = word(comfortFoamDict.lookup("cellSet"));
    }
    if (dictSetFieldsFound)
    {
        const word setFieldsName(word(comfortFoamDict.lookup("setFields")));
        if (dictCellSetFound && setFieldsName != dictCellSetName)
        {
            FatalErrorInFunction
                << "Both dictionary entries cellSet and setFields are defined with "
                << "different names: " << dictCellSetName << " vs " << setFieldsName << nl
                << "Use only one of these entries (or use the same name)."
                << abort(FatalError);
        }
        dictCellSetName = setFieldsName;
    }

    const bool cliAnySelection =
        cliCellZoneFound || !cliCellSetName.empty();
    const bool dictAnySelection =
        dictCellZoneFound || !dictCellSetName.empty();

    if (cliAnySelection && dictAnySelection)
    {
        Info<< "Both CLI and comfortFoamDict region selections are specified. "
            << "Using CLI selection." << endl;
    }

    // CLI selection has priority over dictionary selection
    bool useCellSet = false;
    bool useCellZone = false;
    word cellSetName;
    word cellZoneName;

    if (cliAnySelection)
    {
        useCellSet = !cliCellSetName.empty();
        useCellZone = cliCellZoneFound;
        cellSetName = cliCellSetName;
        if (cliCellZoneFound)
        {
            cellZoneName = args.get<word>("cellZone");
        }
    }
    else
    {
        useCellSet = !dictCellSetName.empty();
        useCellZone = dictCellZoneFound;
        cellSetName = dictCellSetName;
        if (dictCellZoneFound)
        {
            cellZoneName = word(comfortFoamDict.lookup("cellZone"));
        }
    }

    if (useCellSet && useCellZone)
    {
        FatalErrorInFunction
            << "Both cellSet and cellZone region selections are specified. "
            << "Use only one."
            << abort(FatalError);
    }

    // Check for cellSet option
    if (useCellSet)
    {
        Info<< "Loading cellSet: " << cellSetName << endl;
        
        cellSet selectedCells(mesh, cellSetName);

        const labelList rawCells = selectedCells.toc();
        cellsToAnalyze =
            normalizeSelectedCells(mesh, rawCells, "cellSet", cellSetName);

        const label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        if (nCells == 0)
        {
            FatalErrorInFunction
                << "cellSet " << cellSetName
                << " is empty, invalid, or not available on this case."
                << abort(FatalError);
        }

        Info<< "Analyzing " << nCells 
            << " cells from cellSet " << cellSetName << endl;
        selectionType = AnalysisRegionType::CellSet;
        selectionName = cellSetName;
    }
    // Check for cellZone option
    else if (useCellZone)
    {
        Info<< "Loading cellZone: " << cellZoneName << endl;
        
        const cellZoneMesh& cellZones = mesh.cellZones();
        label zoneID = cellZones.findZoneID(cellZoneName);
        
        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneName << " does not exist" << nl
                << "Available cellZones: " << cellZones.names()
                << abort(FatalError);
        }
        
        const cellZone& cz = cellZones[zoneID];
        const labelList rawCells = cz;
        cellsToAnalyze =
            normalizeSelectedCells(mesh, rawCells, "cellZone", cellZoneName);

        const label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        if (nCells == 0)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneName
                << " is empty after parallel label mapping."
                << abort(FatalError);
        }

        Info<< "Analyzing " << nCells 
            << " cells from cellZone " << cellZoneName << endl;
        selectionType = AnalysisRegionType::CellZone;
        selectionName = cellZoneName;
    }
    // Use entire mesh
    else
    {
        cellsToAnalyze.setSize(mesh.nCells());
        forAll(cellsToAnalyze, i)
        {
            cellsToAnalyze[i] = i;
        }
        
        label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        Info<< "Analyzing entire mesh: " << nCells << " cells" << endl;
        selectionType = AnalysisRegionType::EntireMesh;
        selectionName.clear();
    }
    
    return cellsToAnalyze;
}

struct AnnexDValidationCase
{
    scalar airTemperature;
    scalar radiantTemperature;
    scalar relativeAirVelocity;
    scalar relativeHumidity;
    scalar met;
    scalar clo;
    scalar expectedPMV;
    scalar expectedPPD;
};

Foam::label runISO7730Validation()
{
    static const AnnexDValidationCase cases[] =
    {
        {22.0, 22.0, 0.10, 60.0, 1.2, 0.5, -0.75374, 16.965},
        {27.0, 27.0, 0.10, 60.0, 1.2, 0.5,  0.76523, 17.336},
        {27.0, 27.0, 0.30, 60.0, 1.2, 0.5,  0.43456,  8.939},
        {23.5, 25.5, 0.10, 60.0, 1.2, 0.5, -0.01499,  5.005},
        {23.5, 25.5, 0.30, 60.0, 1.2, 0.5, -0.55431, 11.433},
        {19.0, 19.0, 0.10, 40.0, 1.2, 1.0, -0.60129, 12.581},
        {23.5, 23.5, 0.10, 60.0, 1.2, 1.0,  0.48898,  9.996},
        {23.5, 23.5, 0.30, 40.0, 1.2, 1.0,  0.11843,  5.291},
        {23.0, 21.0, 0.10, 40.0, 1.2, 1.0,  0.05125,  5.054},
        {23.0, 21.0, 0.30, 40.0, 1.2, 1.0, -0.16683,  5.577},
        {22.0, 22.0, 0.10, 60.0, 1.6, 0.5,  0.04591,  5.044},
        {27.0, 27.0, 0.10, 60.0, 1.6, 0.5,  1.17102, 33.843},
        {27.0, 27.0, 0.30, 60.0, 1.6, 0.5,  0.95095, 24.101}
    };

    const scalar pmvTolerance = 1e-3;
    const scalar ppdTolerance = 0.05;
    label failures = 0;

    Info<< nl
        << "ISO 7730:2025 Annex D regression validation" << nl
        << "------------------------------------------------" << endl;

    for (label caseI = 0; caseI < 13; ++caseI)
    {
        const AnnexDValidationCase& test = cases[caseI];
        const scalar vapourPressure = calculateWaterVapourPressure
        (
            test.airTemperature,
            test.relativeHumidity
        );
        bool converged = false;
        const scalar pmv = calculatePMV
        (
            test.airTemperature,
            test.relativeAirVelocity,
            vapourPressure,
            test.radiantTemperature,
            test.met*ComfortConstants::baseMetabolicRate,
            test.clo,
            0.0,
            converged
        );
        const scalar ppd = calculatePPD(pmv);
        const bool passed =
            converged
         && mag(pmv - test.expectedPMV) <= pmvTolerance
         && mag(ppd - test.expectedPPD) <= ppdTolerance;

        Info<< "Case " << caseI + 1 << ": PMV=" << pmv
            << ", PPD=" << ppd << " % -- "
            << (passed ? "PASS" : "FAIL") << endl;

        if (!passed)
        {
            ++failures;
            Info<< "  expected PMV=" << test.expectedPMV
                << ", PPD=" << test.expectedPPD << " %" << endl;
        }
    }

    const scalar fieldTemperatureC =
        295.15 - ComfortConstants::physicalKelvinOffset;
    const scalar fieldOperativeTemperature =
        calculateOperativeTemperature(fieldTemperatureC, 22.0, 0.1);
    const bool physicalKelvinPassed =
        mag(fieldTemperatureC - 22.0) < 1e-12
     && mag(fieldOperativeTemperature - 295.15) < 1e-12;
    Info<< "Physical-kelvin field conversion: "
        << (physicalKelvinPassed ? "PASS" : "FAIL") << endl;
    if (!physicalKelvinPassed)
    {
        ++failures;
    }

    const scalar initialTcl = calculateInitialClothingSurfaceTemperature
    (
        295.0,
        22.0,
        0.155*0.5
    );
    const bool initialTclPassed = mag(initialTcl - 301.4299109933617) < 1e-9;
    Info<< "Annex D line 250: "
        << (initialTclPassed ? "PASS" : "FAIL") << endl;
    if (!initialTclPassed)
    {
        ++failures;
    }

    const scalar dr = calculateDraughtRating(22.0, 0.1, 0.4);
    const bool drPassed = mag(dr - 8.653357061761142) < 1e-9;
    Info<< "Draught-rate turbulence-percent conversion: "
        << (drPassed ? "PASS" : "FAIL") << endl;
    if (!drPassed)
    {
        ++failures;
    }

    bool clothingCorrectionApplicable = false;
    const scalar corrected03 = calculateResultantClothingInsulation
    (
        0.3,
        0.5,
        0.2,
        clothingCorrectionApplicable
    );
    const bool dynamicClothingPassed =
        clothingCorrectionApplicable
     && mag(corrected03 - 0.27354745) < 1e-6;
    Info<< "ISO 9920 dynamic-clothing correction: "
        << (dynamicClothingPassed ? "PASS" : "FAIL") << endl;
    if (!dynamicClothingPassed)
    {
        ++failures;
    }

    const scalar correctedLowerVelocity =
        calculateResultantClothingInsulation
        (
            1.0,
            0.0,
            0.2,
            clothingCorrectionApplicable
        );
    const bool clothingVelocityClampPassed =
        clothingCorrectionApplicable
     && mag(correctedLowerVelocity - 0.9086965858) < 1e-6;
    Info<< "Dynamic-clothing 0.15 m/s lower limit: "
        << (clothingVelocityClampPassed ? "PASS" : "FAIL") << endl;
    if (!clothingVelocityClampPassed)
    {
        ++failures;
    }

    const bool categoryBoundsPassed =
        analyzeWholeBodyCategory(0.2, calculatePPD(0.2))
            != ComfortCategory::I
     && analyzeWholeBodyCategory(-0.2, calculatePPD(-0.2))
            != ComfortCategory::I;
    Info<< "Strict Annex A category bounds: "
        << (categoryBoundsPassed ? "PASS" : "FAIL") << endl;
    if (!categoryBoundsPassed)
    {
        ++failures;
    }

    const scalar verticalPD =
        calculateVerticalTemperatureDissatisfaction(1.0);
    const scalar floorPD = calculateFloorTemperatureDissatisfaction(24.0);
    scalar k1 = 0.0;
    scalar k2 = 0.0;
    scalar k3 = 0.0;
    scalar asymmetryLimit = 0.0;
    const bool radiantConstantsValid = radiantAsymmetryConstants
    (
        "warmCeiling",
        k1,
        k2,
        k3,
        asymmetryLimit
    );
    const scalar radiantPD = calculateRadiantAsymmetryDissatisfaction
    (
        2.0,
        k1,
        k2,
        k3
    );
    const bool localModelsPassed =
        radiantConstantsValid
     && asymmetryLimit == 23.0
     && mag(verticalPD - 0.736225) < 1e-5
     && mag(floorPD - 5.52882) < 1e-5
     && mag(radiantPD - 1.36253) < 1e-5
     && analyzeOverallComfortCategory
        (
            0.0,
            5.0,
            8.0,
            verticalPD,
            floorPD,
            radiantPD,
            true
        ) == ComfortCategory::I;
    Info<< "Local-discomfort models and Category I: "
        << (localModelsPassed ? "PASS" : "FAIL") << endl;
    if (!localModelsPassed)
    {
        ++failures;
    }

    Info<< "------------------------------------------------" << nl
        << (failures == 0 ? "VALIDATION PASSED" : "VALIDATION FAILED")
        << " (" << failures << " failure(s))" << nl << endl;

    return failures == 0 ? 0 : 1;
}

// * * * * * * * * * * * * * * * * Main Program * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    
    argList::addBoolOption
    (
        "validate",
        "Run all ISO 7730:2025 Annex D regression tests"
    );
    
    argList::addOption
    (
        "cellSet",
        "name",
        "Analyze only cells in the specified cellSet instead of entire mesh"
    );

    argList::addOption
    (
        "setFields",
        "name",
        "Alias for -cellSet (analyze only cells in the specified cellSet)"
    );
    
    argList::addOption
    (
        "cellZone", 
        "name",
        "Analyze only cells in the specified cellZone instead of entire mesh"
    );
    
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    
    if (args.found("validate"))
    {
        return runISO7730Validation();
    }

    // Get times list
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    if (timeDirs.size() > 1)
    {
        WarningInFunction
            << "Multiple time directories were selected. comfortFoam evaluates "
            << "each snapshot independently; it does not create the one-hour "
            << "time-weighted inputs required by ISO 7730 for lightly varying "
            << "conditions." << endl;
    }

    #include "createNamedMesh.H"

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        // Load temperature field
        volScalarField T(THeader, mesh);
        const fvPatchList& patches = T.mesh().boundary();

        validateInputParameters(met, clo, wme, RH1);

        if (applyDynamicClothing && clothingAlreadyResultant)
        {
            FatalErrorInFunction
                << "applyDynamicClothing and clothingAlreadyResultant cannot "
                << "both be true. Choose exactly one clothing-input mode."
                << abort(FatalError);
        }

        const bool clothingModeDeclared =
            applyDynamicClothing || clothingAlreadyResultant;
        if (!clothingModeDeclared)
        {
            WarningInFunction
                << "No ISO 7730:2025 clothing-input mode is declared. "
                << "The legacy clo value will be used, but ISO7730Valid "
                << "will be zero. Set applyDynamicClothing or "
                << "clothingAlreadyResultant to true." << endl;
        }

        if
        (
            configuredWalkingSpeed > 1.2
         ||
            (
                configuredWalkingSpeed < 0.0
             && configuredWalkingSpeed != -1.0
            )
        )
        {
            FatalErrorInFunction
                << "walkingSpeed must be -1 (automatic) or between 0 and "
                << "1.2 m/s: " << configuredWalkingSpeed
                << abort(FatalError);
        }

        scalar radiantK1 = 0.0;
        scalar radiantK2 = 0.0;
        scalar radiantK3 = 0.0;
        scalar radiantAsymmetryLimit = 0.0;
        if
        (
            !radiantAsymmetryConstants
            (
                radiantAsymmetryType,
                radiantK1,
                radiantK2,
                radiantK3,
                radiantAsymmetryLimit
            )
        )
        {
            FatalErrorInFunction
                << "Invalid radiantAsymmetryType '" << radiantAsymmetryType
                << "'. Valid values are warmCeiling, coolWall, coolCeiling "
                << "and warmWall."
                << abort(FatalError);
        }

        // Check for humidity fields from supported solver conventions.
        bool humidityAvailable = false;
        autoPtr<volScalarField> humidityField;
        word humiditySource("comfortFoamDict RH");

        const wordList humidityCandidates
        ({
            "thermo:relHum",
            "thermoRelHum",
            "relHum",
            "RH",
            "relativeHumidity"
        });

        forAll(humidityCandidates, i)
        {
            IOobject humidityHeader
            (
                humidityCandidates[i],
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT
            );

            if (humidityHeader.typeHeaderOk<volScalarField>())
            {
                humidityField.reset(new volScalarField(humidityHeader, mesh));
                humidityAvailable = true;
                humiditySource = humidityCandidates[i];
                break;
            }
        }

        autoPtr<volScalarField> waterVapourPressureField;
        if (waterVapourPressureHeader.typeHeaderOk<volScalarField>())
        {
            waterVapourPressureField.reset
            (
                new volScalarField(waterVapourPressureHeader, mesh)
            );
            Info<< "Using water-vapour pressure field: waterVapourPressure"
                << endl;
        }
        else if (configuredWaterVapourPressure >= 0.0)
        {
            Info<< "Using constant water-vapour pressure: "
                << configuredWaterVapourPressure << " Pa" << endl;
        }

        if
        (
            !waterVapourPressureField.valid()
         && configuredWaterVapourPressure < 0.0
         && humidityAvailable
        )
        {
            Info<< "Using humidity field: " << humiditySource << endl;
        }
        else if
        (
            !waterVapourPressureField.valid()
         && configuredWaterVapourPressure < 0.0
        )
        {
            Info<< "No humidity field found - using default value: " << RH1 << "%" << endl;
        }

        // Turbulence intensity requires only k.
        bool turbulenceAvailable = false;
        autoPtr<volScalarField> kField;

        if (kHeader.typeHeaderOk<volScalarField>())
        {
            kField.reset(new volScalarField(kHeader, mesh));
            turbulenceAvailable = true;
            
            if (epsilonHeader.typeHeaderOk<volScalarField>())
            {
                Info<< "Turbulence source: k (k-epsilon model detected)"
                    << endl;
            }
            else if (omegaHeader.typeHeaderOk<volScalarField>())
            {
                Info<< "Turbulence source: k (k-omega model detected)"
                    << endl;
            }
            else
            {
                Info<< "Turbulence source: k" << endl;
            }
        }
        else
        {
            Info<< "No turbulence fields available - assuming Tu = 40%" << endl;
        }

        // Load velocity field
        volVectorField U(UHeader, mesh);

        autoPtr<volScalarField> mrtField;
        if (MRTHeader.typeHeaderOk<volScalarField>())
        {
            mrtField.reset(new volScalarField(MRTHeader, mesh));
            Info<< "Using mean radiant temperature field: MRT" << endl;
        }

        autoPtr<volScalarField> relativeAirVelocityField;
        if (relativeAirVelocityHeader.typeHeaderOk<volScalarField>())
        {
            relativeAirVelocityField.reset
            (
                new volScalarField(relativeAirVelocityHeader, mesh)
            );
            Info<< "Using relative-air-velocity field: relativeAirVelocity"
                << endl;
        }
        else if (configuredRelativeAirVelocity >= 0.0)
        {
            Info<< "Using constant relative air velocity: "
                << configuredRelativeAirVelocity << " m/s" << endl;
        }

        autoPtr<volScalarField> verticalTemperatureDifferenceField;
        if
        (
            verticalTemperatureDifferenceHeader
               .typeHeaderOk<volScalarField>()
        )
        {
            verticalTemperatureDifferenceField.reset
            (
                new volScalarField(verticalTemperatureDifferenceHeader, mesh)
            );
        }

        autoPtr<volScalarField> floorTemperatureField;
        if (floorTemperatureHeader.typeHeaderOk<volScalarField>())
        {
            floorTemperatureField.reset
            (
                new volScalarField(floorTemperatureHeader, mesh)
            );
        }

        autoPtr<volScalarField> radiantTemperatureAsymmetryField;
        if
        (
            radiantTemperatureAsymmetryHeader
               .typeHeaderOk<volScalarField>()
        )
        {
            radiantTemperatureAsymmetryField.reset
            (
                new volScalarField(radiantTemperatureAsymmetryHeader, mesh)
            );
        }
        
        // Get cells to analyze (cellSet/cellZone from CLI or comfortFoamDict)
        AnalysisRegionType selectionType;
        word selectionName;
        labelList cellsToAnalyze =
            getCellsToAnalyze(mesh, args, comfortFoamDict, selectionType, selectionName);
        
        // Calculate averages for analysis region
        scalar avgVelocityMag, avgTemperature, totalAnalysisVolume;
        calculateVolumeWeightedAverages(mesh, cellsToAnalyze, avgVelocityMag, avgTemperature, totalAnalysisVolume);

        // Pre-calculate constants
        const scalar metRate = met * ComfortConstants::baseMetabolicRate;
        const scalar wmeRate = wme * ComfortConstants::baseMetabolicRate;
        const scalar walkingSpeed = applyDynamicClothing
            ? calculateWalkingSpeed(metRate, configuredWalkingSpeed)
            : Foam::max(0.0, configuredWalkingSpeed);

        if (applyDynamicClothing && (clo < 0.0 || clo >= 1.4))
        {
            WarningInFunction
                << "Dynamic clothing correction is defined for 0 <= clo < "
                << "1.4. The configured value is " << clo << " clo. "
                << "Use clothingAlreadyResultant with a pre-corrected value "
                << "for this case." << endl;
        }

        // Volume-weighted averages for output
        scalar volumeWeightedPMV(0);
        scalar volumeWeightedPPD(0);
        scalar volumeWeightedDR(0);
        scalar volumeWeightedTOp(0);
        scalar volumeWeightedRadTemp(0);
        scalar volumeWeightedRH(0);
        scalar volumeWeightedTu(0);
        scalar volumeWeightedRelativeVelocity(0);
        scalar volumeWeightedResultantClo(0);
        scalar totalVolume(0);

        label invalidPMVCells = 0;
        label invalidDRCells = 0;
        label invalidLocalCells = 0;
        label missingLocalCells = 0;
        label nonConvergedCells = 0;
        bool anyWholeBodyNoCategory = false;
        bool anyOverallNoCategory = false;
        bool anyOverallIncomplete = false;
        label worstWholeBodyCategory = 1;
        label worstOverallCategory = 1;

        // Radiation source selection must be rank-consistent in parallel.
        const bool useMRTField = mrtField.valid();
        const bool useGField = !useMRTField && G.headerOk();
        const bool useQrField =
            !useMRTField && !useGField
         && qrHeader.typeHeaderOk<volScalarField>();
        const bool useIDefaultField =
            !useMRTField && !useGField && !useQrField
         && IDefaultHeader.typeHeaderOk<volScalarField>();

        autoPtr<volScalarField> qrField;
        autoPtr<volScalarField> IDefaultField;
        scalar fallbackRadTemp = 0.0;

        if (useQrField)
        {
            qrField.reset(new volScalarField(qrHeader, mesh));
        }
        else if (useIDefaultField)
        {
            IDefaultField.reset(new volScalarField(IDefaultHeader, mesh));
        }
        else if (!useMRTField && !useGField)
        {
            // Uses gSum internally, so evaluate once per time-step (not per cell).
            fallbackRadTemp = calculateRadiationTemperature(mesh, patches);
        }

        // Main calculation loop - only over selected cells
        forAll(cellsToAnalyze, i)
        {
            const label cellI = cellsToAnalyze[i];
            const scalar cellTemp = T[cellI];
            const scalar cellTempC =
                cellTemp - ComfortConstants::physicalKelvinOffset;
            const scalar localAirVelocity = mag(U[cellI]);

            // Determine radiation temperature from available radiation fields
            scalar cellRadTemp;
            if (useMRTField)
            {
                cellRadTemp =
                    mrtField()[cellI]
                  - ComfortConstants::physicalKelvinOffset;
            }
            else if (useGField)
            {
                scalar gValue = Foam::max(0.0, Foam::min(50000.0, G[cellI]));
                cellRadTemp = Foam::pow
                (
                    gValue
                   /(4.0*ComfortConstants::stefanBoltzmannConstant),
                    0.25
                ) - ComfortConstants::physicalKelvinOffset;
            }
            else if (qrField.valid())
            {
                const scalar emissivity = 0.9;
                const scalar qrValue = qrField()[cellI];
                const scalar T4_rad = Foam::pow(cellTemp, 4)
                    - qrValue
                     /(emissivity*ComfortConstants::stefanBoltzmannConstant);
                if (T4_rad > 0)
                {
                    cellRadTemp =
                        Foam::pow(T4_rad, 0.25)
                      - ComfortConstants::physicalKelvinOffset;
                }
                else
                {
                    cellRadTemp = cellTempC;
                }
            }
            else if (IDefaultField.valid())
            {
                scalar gValue = 4.0*constant::mathematical::pi
                    *IDefaultField()[cellI];
                gValue = Foam::max(0.0, Foam::min(50000.0, gValue));
                cellRadTemp = Foam::pow
                (
                    gValue
                   /(4.0*ComfortConstants::stefanBoltzmannConstant),
                    0.25
                ) - ComfortConstants::physicalKelvinOffset;
            }
            else
            {
                cellRadTemp = fallbackRadTemp;
            }

            scalar cellRH = RH1;
            scalar cellVapourPressure = -1.0;
            if (waterVapourPressureField.valid())
            {
                cellVapourPressure = waterVapourPressureField()[cellI];
                cellRH = calculateRelativeHumidity
                (
                    cellTempC,
                    cellVapourPressure
                );
            }
            else if (configuredWaterVapourPressure >= 0.0)
            {
                cellVapourPressure = configuredWaterVapourPressure;
                cellRH = calculateRelativeHumidity
                (
                    cellTempC,
                    cellVapourPressure
                );
            }
            else if (humidityAvailable)
            {
                const scalar rawHumidity = humidityField()[cellI];
                const bool knownFractionSource =
                    humiditySource == "thermo:relHum"
                 || humiditySource == "thermoRelHum";
                const bool knownPercentSource =
                    humiditySource == "RH"
                 || humiditySource == "relativeHumidity";

                if
                (
                    knownFractionSource
                 || (!knownPercentSource && rawHumidity <= 1.0)
                )
                {
                    cellRH = rawHumidity*100.0;
                }
                else
                {
                    cellRH = rawHumidity;
                }
            }

            cellVapourPressure = cellVapourPressure >= 0.0
                ? cellVapourPressure
                : calculateWaterVapourPressure(cellTempC, cellRH);

            scalar relativeAirVelocity = localAirVelocity + walkingSpeed;
            if (relativeAirVelocityField.valid())
            {
                relativeAirVelocity = relativeAirVelocityField()[cellI];
            }
            else if (configuredRelativeAirVelocity >= 0.0)
            {
                relativeAirVelocity = configuredRelativeAirVelocity;
            }
            const scalar modelRelativeAirVelocity =
                Foam::max(relativeAirVelocity, 0.0);

            bool clothingCorrectionApplicable = true;
            scalar resultantClo = clo;
            if (applyDynamicClothing)
            {
                resultantClo = calculateResultantClothingInsulation
                (
                    clo,
                    modelRelativeAirVelocity,
                    walkingSpeed,
                    clothingCorrectionApplicable
                );
            }
            const bool clothingInputCompliant =
                clothingModeDeclared && clothingCorrectionApplicable;

            const scalar turbulenceIntensity =
                calculateTurbulentIntensityFraction
            (
                localAirVelocity,
                turbulenceAvailable && kField.valid() ? kField()[cellI] : 0.0,
                turbulenceAvailable && kField.valid()
            );

            bool clothingIterationConverged = false;
            const scalar pmv = calculatePMV
            (
                cellTempC,
                modelRelativeAirVelocity,
                cellVapourPressure,
                cellRadTemp,
                metRate,
                resultantClo,
                wmeRate,
                clothingIterationConverged
            );
            const scalar ppd = calculatePPD(pmv);
            const scalar dr = calculateDraughtRating
            (
                cellTempC,
                localAirVelocity,
                turbulenceIntensity
            );
            const scalar tOp = calculateOperativeTemperature
            (
                cellTempC,
                cellRadTemp,
                localAirVelocity
            );

            const bool pmvApplicable = isPMVModelApplicable
            (
                cellTempC,
                cellRadTemp,
                relativeAirVelocity,
                cellVapourPressure,
                resultantClo,
                pmv,
                clothingIterationConverged,
                clothingInputCompliant
            );
            const bool drApplicable = isDraughtModelApplicable
            (
                cellTempC,
                localAirVelocity,
                turbulenceIntensity,
                met
            );

            const bool verticalAvailable =
                verticalTemperatureDifferenceField.valid()
             || configuredVerticalTemperatureDifference >= 0.0;
            const scalar verticalDifference =
                verticalTemperatureDifferenceField.valid()
              ? verticalTemperatureDifferenceField()[cellI]
              : configuredVerticalTemperatureDifference;
            const bool verticalValid =
                verticalAvailable
             && verticalDifference >= 0.0
             && verticalDifference < 8.0;
            const scalar verticalPD = verticalValid
                ? calculateVerticalTemperatureDissatisfaction
                  (
                      verticalDifference
                  )
                : -1.0;

            const bool floorAvailable =
                floorTemperatureField.valid()
             || configuredFloorTemperature >= 0.0;
            const scalar floorTemperatureC = floorTemperatureField.valid()
                ? floorTemperatureField()[cellI]
                  - ComfortConstants::physicalKelvinOffset
                : configuredFloorTemperature;
            const bool floorValid = floorAvailable;
            const scalar floorPD = floorValid
                ? calculateFloorTemperatureDissatisfaction(floorTemperatureC)
                : -1.0;

            const bool radiantAvailable =
                radiantTemperatureAsymmetryField.valid()
             || configuredRadiantTemperatureAsymmetry >= 0.0;
            const scalar radiantAsymmetry =
                radiantTemperatureAsymmetryField.valid()
              ? radiantTemperatureAsymmetryField()[cellI]
              : configuredRadiantTemperatureAsymmetry;
            const bool radiantValid =
                radiantAvailable
             && radiantAsymmetry >= 0.0
             && radiantAsymmetry < radiantAsymmetryLimit;
            const scalar radiantPD = radiantValid
                ? calculateRadiantAsymmetryDissatisfaction
                  (
                      radiantAsymmetry,
                      radiantK1,
                      radiantK2,
                      radiantK3
                  )
                : -1.0;

            const bool localInputsAvailable =
                verticalAvailable && floorAvailable && radiantAvailable;
            const bool localInputsValid =
                verticalValid && floorValid && radiantValid;
            const bool completeLocalAssessment =
                drApplicable && localInputsAvailable && localInputsValid;

            const ComfortCategory wholeBodyCategory = pmvApplicable
                ? analyzeWholeBodyCategory(pmv, ppd)
                : ComfortCategory::None;
            const ComfortCategory overallCategory = !pmvApplicable
                ? ComfortCategory::None
                :
                  (
                      completeLocalAssessment
                    ? analyzeOverallComfortCategory
                      (
                          pmv,
                          ppd,
                          dr,
                          verticalPD,
                          floorPD,
                          radiantPD,
                          true
                      )
                    :
                      (
                          wholeBodyCategory == ComfortCategory::IV
                        ? ComfortCategory::IV
                        : ComfortCategory::Incomplete
                      )
                  );

            PMV[cellI] = pmv;
            PPD[cellI] = ppd;
            DR[cellI] = dr;
            TOp[cellI] = tOp;
            PDVertical[cellI] = verticalPD;
            PDFloor[cellI] = floorPD;
            PDRadiantAsymmetry[cellI] = radiantPD;
            ISO7730Valid[cellI] = pmvApplicable ? 1.0 : 0.0;
            ISO7730WholeBodyCategory[cellI] =
                comfortCategoryValue(wholeBodyCategory);
            ISO7730Category[cellI] =
                comfortCategoryValue(overallCategory);

            if (!pmvApplicable)
            {
                ++invalidPMVCells;
            }
            if (!drApplicable)
            {
                ++invalidDRCells;
            }
            if (!clothingIterationConverged)
            {
                ++nonConvergedCells;
            }
            if (!localInputsAvailable)
            {
                ++missingLocalCells;
            }
            else if (!localInputsValid)
            {
                ++invalidLocalCells;
            }

            if (wholeBodyCategory == ComfortCategory::None)
            {
                anyWholeBodyNoCategory = true;
            }
            else
            {
                worstWholeBodyCategory = Foam::max
                (
                    worstWholeBodyCategory,
                    static_cast<label>(wholeBodyCategory)
                );
            }

            if (overallCategory == ComfortCategory::Incomplete)
            {
                anyOverallIncomplete = true;
            }
            else if (overallCategory == ComfortCategory::None)
            {
                // An inapplicable PMV is indeterminate, not a definite
                // category failure. Track only evaluated no-category cells.
                if (pmvApplicable)
                {
                    anyOverallNoCategory = true;
                }
            }
            else
            {
                worstOverallCategory = Foam::max
                (
                    worstOverallCategory,
                    static_cast<label>(overallCategory)
                );
            }

            // Accumulate volume-weighted averages
            const scalar cellVolume = mesh.V()[cellI];
            volumeWeightedPMV += pmv*cellVolume;
            volumeWeightedPPD += ppd*cellVolume;
            volumeWeightedDR += dr*cellVolume;
            volumeWeightedTOp += tOp*cellVolume;
            volumeWeightedRadTemp += cellRadTemp*cellVolume;
            volumeWeightedRH += cellRH*cellVolume;
            volumeWeightedTu += turbulenceIntensity*cellVolume;
            volumeWeightedRelativeVelocity += relativeAirVelocity*cellVolume;
            volumeWeightedResultantClo += resultantClo*cellVolume;
            totalVolume += cellVolume;
        } // End of forAll(cellsToAnalyze, i)

        // Write fields
        DR.write();
        PMV.write();
        PPD.write();
        TOp.write();
        PDVertical.write();
        PDFloor.write();
        PDRadiantAsymmetry.write();
        ISO7730Valid.write();
        ISO7730WholeBodyCategory.write();
        ISO7730Category.write();

        // Reduce parallel values
        reduce(volumeWeightedPMV, sumOp<scalar>());
        reduce(volumeWeightedPPD, sumOp<scalar>());
        reduce(volumeWeightedDR, sumOp<scalar>());
        reduce(volumeWeightedTOp, sumOp<scalar>());
        reduce(volumeWeightedRadTemp, sumOp<scalar>());
        reduce(volumeWeightedRH, sumOp<scalar>());
        reduce(volumeWeightedTu, sumOp<scalar>());
        reduce(volumeWeightedRelativeVelocity, sumOp<scalar>());
        reduce(volumeWeightedResultantClo, sumOp<scalar>());
        reduce(totalVolume, sumOp<scalar>());
        reduce(invalidPMVCells, sumOp<label>());
        reduce(invalidDRCells, sumOp<label>());
        reduce(invalidLocalCells, sumOp<label>());
        reduce(missingLocalCells, sumOp<label>());
        reduce(nonConvergedCells, sumOp<label>());
        reduce(anyWholeBodyNoCategory, orOp<bool>());
        reduce(anyOverallNoCategory, orOp<bool>());
        reduce(anyOverallIncomplete, orOp<bool>());
        reduce(worstWholeBodyCategory, maxOp<label>());
        reduce(worstOverallCategory, maxOp<label>());

        if (totalVolume <= SMALL)
        {
            FatalErrorInFunction
                << "The selected analysis region has zero volume."
                << abort(FatalError);
        }
        
        // Calculate final averages
        scalar avgPMV = volumeWeightedPMV / totalVolume;
        scalar avgPPD = volumeWeightedPPD / totalVolume;
        scalar avgDR = volumeWeightedDR / totalVolume;
        scalar avgTOp = volumeWeightedTOp / totalVolume;
        scalar avgRadTemp = volumeWeightedRadTemp / totalVolume;
        scalar avgRH = volumeWeightedRH / totalVolume;
        scalar avgTu = volumeWeightedTu / totalVolume;
        scalar avgRelativeVelocity =
            volumeWeightedRelativeVelocity/totalVolume;
        scalar avgResultantClo = volumeWeightedResultantClo/totalVolume;

        const ComfortCategory worstWholeBody = anyWholeBodyNoCategory
            ? ComfortCategory::None
            : static_cast<ComfortCategory>(worstWholeBodyCategory);

        const char* overallCategoryText = "No category";
        if (anyOverallNoCategory)
        {
            overallCategoryText = "No category";
        }
        else if (invalidPMVCells > 0)
        {
            overallCategoryText =
                "Not determined (one or more inputs are outside "
                "ISO 7730 applicability)";
        }
        else if (anyOverallIncomplete)
        {
            overallCategoryText =
                "Not determined (local inputs are missing or outside "
                "applicability)";
        }
        else
        {
            overallCategoryText = comfortCategoryToString
            (
                static_cast<ComfortCategory>(worstOverallCategory)
            );
        }

        // Output results
        Info<< nl << "============ THERMAL COMFORT ANALYSIS RESULTS ============" << nl;
        
        // Show analysis region info
        if (selectionType == AnalysisRegionType::CellSet)
        {
            Info<< "Analysis region: cellSet " << selectionName << nl;
        }
        else if (selectionType == AnalysisRegionType::CellZone)
        {
            Info<< "Analysis region: cellZone " << selectionName << nl;
        }
        else
        {
            Info<< "Analysis region: entire mesh" << nl;
        }
        
        label totalCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        Info<< "Analyzed cells: " << totalCells << nl
            << "Analysis volume: " << totalVolume << " m^3" << nl
            << nl
            << "Mean radiation temperature:     " << avgRadTemp << " degrees C" << nl
            << "Average air temperature:        "
            << avgTemperature - ComfortConstants::physicalKelvinOffset
            << " degrees C" << nl
            << "Average air velocity:           " << avgVelocityMag << " m/s" << nl
            << "Average relative air velocity:  " << avgRelativeVelocity << " m/s" << nl
            << "Average relative humidity:      " << avgRH << " %" << nl
            << "Average resultant clothing:     " << avgResultantClo << " clo" << nl
            << "Average operative temperature:  "
            << avgTOp - ComfortConstants::physicalKelvinOffset
            << " degrees C" << nl
            << "Average turbulent intensity:    " << avgTu * 100.0 << " %" << nl
            << nl
            << "COMFORT INDICES:" << nl
            << "PMV (Predicted Mean Vote):      " << avgPMV << nl
            << "PPD (Predicted % Dissatisfied): " << avgPPD << " %" << nl
            << "DR (Draught Rating):            " << avgDR << " %" << nl
            << nl
            << "ISO 7730:2025 APPLICABILITY:" << nl
            << "PMV cells outside applicability: " << invalidPMVCells
            << " / " << totalCells << nl
            << "DR cells outside applicability:  " << invalidDRCells
            << " / " << totalCells << nl
            << "Local input missing in cells:    " << missingLocalCells
            << " / " << totalCells << nl
            << "Local input invalid in cells:    " << invalidLocalCells
            << " / " << totalCells << nl
            << "Clothing iterations not converged: " << nonConvergedCells
            << nl
            << nl
            << "Worst whole-body category: "
            <<
               (
                   invalidPMVCells > 0
                 ? "Not determined (PMV outside applicability)"
                 : comfortCategoryToString(worstWholeBody)
               )
            << nl
            << "Worst overall category:    " << overallCategoryText << nl
            << "=========================================================" << endl;
    
    } // End of forAll(timeDirs, timei)

    return 0;
}

// ************************************************************************* //
