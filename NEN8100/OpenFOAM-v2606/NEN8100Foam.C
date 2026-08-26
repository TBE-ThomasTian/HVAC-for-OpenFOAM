/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Adjusted by FlowSIM Pro, 2026.

License
    This file is part of HVAC-for-OpenFOAM.

    HVAC-for-OpenFOAM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Application
    NEN8100Foam

Description
    Multi-case pedestrian wind-comfort post-processing according to the
    public NEN 8100:2006 comfort and danger criteria.

    The utility calculates volume fields of local velocity amplification for
    matching directional CFD cases and combines them with directional wind
    statistics. It also writes a pedestrian-height mask for ParaView.

    The initial climate reader supports directional Weibull distributions.
    Omitting the climate dictionary produces amplification factors only.

Notes
    - This initial implementation is serial.
    - Directional cases must use matching cell topology and coordinates.
    - A result from this utility is an engineering assessment, not automatic
      certification of the CFD setup or weather data.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sampledPlane.H"
#include "OSspecific.H"
#include "PtrList.H"

#include <cmath>

namespace Foam
{

bool finiteValue(const scalar value)
{
    return std::isfinite(value);
}


bool finiteValue(const vector& value)
{
    return
        std::isfinite(value.x())
     && std::isfinite(value.y())
     && std::isfinite(value.z());
}

struct DirectionSpec
{
    scalar windFrom;
    fileName caseDir;
    scalar Uref;
    scalar referenceHeight;
    word timeSelection;
    word region;

    DirectionSpec()
    :
        windFrom(0),
        caseDir(),
        Uref(0),
        referenceHeight(0),
        timeSelection("latest"),
        region(polyMesh::defaultRegion)
    {}
};


struct ClimateSector
{
    scalar windFrom;
    scalar frequency;
    scalar scale;
    scalar shape;

    ClimateSector()
    :
        windFrom(0),
        frequency(0),
        scale(0),
        shape(0)
    {}
};


struct MeshSignature
{
    label nPoints;
    label nFaces;
    label nInternalFaces;
    vectorField cellCentres;
    scalarField cellVolumes;

    explicit MeshSignature(const fvMesh& mesh)
    :
        nPoints(mesh.nPoints()),
        nFaces(mesh.nFaces()),
        nInternalFaces(mesh.nInternalFaces()),
        cellCentres(mesh.cellCentres()),
        cellVolumes(mesh.V())
    {
        forAll(cellCentres, cellI)
        {
            if
            (
               !finiteValue(cellCentres[cellI])
             || !finiteValue(cellVolumes[cellI])
             || cellVolumes[cellI] <= VSMALL
            )
            {
                FatalErrorInFunction
                    << "Reference mesh contains an invalid cell centre or "
                    << "volume at cell " << cellI
                    << exit(FatalError);
            }
        }
    }
};


scalar normalizeDirection(const scalar angle)
{
    scalar result = std::fmod(angle, 360.0);
    if (result < 0)
    {
        result += 360.0;
    }
    if (mag(result) < SMALL || mag(result - 360.0) < SMALL)
    {
        result = 0;
    }
    return result;
}


scalar directionDifference(const scalar a, const scalar b)
{
    const scalar direct = mag(normalizeDirection(a) - normalizeDirection(b));
    return min(direct, 360.0 - direct);
}


word directionFieldName(const word& prefix, const scalar windFrom)
{
    return
        prefix
      + "Gamma_"
      + word::printf("%.12g", normalizeDirection(windFrom));
}


label comfortClass(const scalar probabilityPercent)
{
    if (probabilityPercent < 2.5) return 1;  // A
    if (probabilityPercent < 5.0) return 2;  // B
    if (probabilityPercent < 10.0) return 3; // C
    if (probabilityPercent < 20.0) return 4; // D
    return 5;                                // E
}


label dangerClass(const scalar probabilityPercent)
{
    // Published NEN tables assign the 0.05--0.30 % interval to limited risk.
    if (probabilityPercent < 0.05) return 0; // no danger
    if (probabilityPercent < 0.30) return 1; // limited risk
    return 2;                                // dangerous
}


label traversingRating(const label cls)
{
    if (cls <= 3) return 0; // good
    if (cls == 4) return 1; // moderate
    return 2;               // poor
}


label strollingRating(const label cls)
{
    if (cls <= 2) return 0;
    if (cls == 3) return 1;
    return 2;
}


label sittingRating(const label cls)
{
    if (cls == 1) return 0;
    if (cls == 2) return 1;
    return 2;
}


scalar weibullExceedance
(
    const scalar localAmplification,
    const scalar localThreshold,
    const scalar scale,
    const scalar shape
)
{
    if (localAmplification <= VSMALL)
    {
        return 0;
    }

    const scalar referenceThreshold = localThreshold/localAmplification;
    const scalar exceedance =
        Foam::exp(-Foam::pow(referenceThreshold/scale, shape));

    if (!finiteValue(exceedance))
    {
        FatalErrorInFunction
            << "Non-finite Weibull exceedance for amplification "
            << localAmplification << ", threshold " << localThreshold
            << ", scale " << scale << " and shape " << shape
            << exit(FatalError);
    }

    return exceedance;
}


fileName resolveCasePath
(
    const fileName& driverCase,
    const fileName& configuredPath
)
{
    fileName result(configuredPath);
    if (!result.isAbsolute())
    {
        result = driverCase/result;
    }
    result.clean();
    return result;
}


void selectCaseTime(Time& caseTime, const word& selection)
{
    const instantList available(caseTime.times());

    if (available.empty())
    {
        FatalErrorInFunction
            << "No time directories found in " << caseTime.globalPath()
            << exit(FatalError);
    }

    label selected = -1;

    if (selection == "latest")
    {
        selected = available.size() - 1;
    }
    else
    {
        forAll(available, timeI)
        {
            if (available[timeI].name() == selection)
            {
                selected = timeI;
                break;
            }
        }
    }

    if (selected < 0)
    {
        FatalErrorInFunction
            << "Requested time '" << selection << "' was not found in "
            << caseTime.globalPath() << nl
            << "Available times: " << available
            << exit(FatalError);
    }

    caseTime.setTime(available[selected], selected);
}


List<DirectionSpec> readDirections
(
    const dictionary& dict,
    const fileName& driverCase
)
{
    PtrList<dictionary> entries(dict.lookup("directions"));

    if (entries.empty())
    {
        FatalIOErrorInFunction(dict)
            << "The directions list is empty"
            << exit(FatalIOError);
    }

    List<DirectionSpec> result(entries.size());

    forAll(entries, directionI)
    {
        const dictionary& entry = entries[directionI];
        DirectionSpec& spec = result[directionI];

        const scalar configuredWindFrom = entry.get<scalar>("windFrom");
        if (!finiteValue(configuredWindFrom))
        {
            FatalIOErrorInFunction(entry)
                << "windFrom must be finite"
                << exit(FatalIOError);
        }

        spec.windFrom = normalizeDirection(configuredWindFrom);
        spec.caseDir = resolveCasePath
        (
            driverCase,
            entry.get<fileName>("case")
        );
        spec.Uref = entry.get<scalar>("Uref");
        spec.referenceHeight = entry.get<scalar>("referenceHeight");
        spec.timeSelection = entry.getOrDefault<word>("time", "latest");
        spec.region = entry.getOrDefault<word>
        (
            "region",
            polyMesh::defaultRegion
        );

        if (!finiteValue(spec.Uref) || spec.Uref <= VSMALL)
        {
            FatalIOErrorInFunction(entry)
                << "Uref must be positive for windFrom " << spec.windFrom
                << exit(FatalIOError);
        }

        if
        (
           !finiteValue(spec.referenceHeight)
         || spec.referenceHeight <= VSMALL
        )
        {
            FatalIOErrorInFunction(entry)
                << "referenceHeight must be positive for windFrom "
                << spec.windFrom << exit(FatalIOError);
        }

        if (!isDir(spec.caseDir))
        {
            FatalIOErrorInFunction(entry)
                << "Directional case directory does not exist: "
                << spec.caseDir << exit(FatalIOError);
        }

        for (label previous = 0; previous < directionI; ++previous)
        {
            if
            (
                directionDifference
                (
                    spec.windFrom,
                    result[previous].windFrom
                ) < SMALL
            )
            {
                FatalIOErrorInFunction(entry)
                    << "Duplicate windFrom direction " << spec.windFrom
                    << exit(FatalIOError);
            }

            if
            (
                directionFieldName("NEN8100", spec.windFrom)
             == directionFieldName
                (
                    "NEN8100",
                    result[previous].windFrom
                )
            )
            {
                FatalIOErrorInFunction(entry)
                    << "Directions " << spec.windFrom << " and "
                    << result[previous].windFrom
                    << " produce the same amplification field name"
                    << exit(FatalIOError);
            }
        }
    }

    return result;
}


List<ClimateSector> readClimate
(
    const dictionary& climateDict,
    const List<DirectionSpec>& directions
)
{
    const word climateType(climateDict.get<word>("type"));
    if (climateType != "directionalWeibull")
    {
        FatalIOErrorInFunction(climateDict)
            << "Unsupported climate type '" << climateType << "'. "
            << "The initial implementation supports directionalWeibull."
            << exit(FatalIOError);
    }

    const scalar referenceHeight =
        climateDict.get<scalar>("referenceHeight");
    const scalar angleTolerance =
        climateDict.getOrDefault<scalar>("directionTolerance", 0.1);
    const scalar calmFraction =
        climateDict.getOrDefault<scalar>("calmFraction", 0);
    const scalar frequencyTolerance =
        climateDict.getOrDefault<scalar>("frequencyTolerance", 1e-8);

    if
    (
       !finiteValue(referenceHeight)
     || !finiteValue(angleTolerance)
     || !finiteValue(calmFraction)
     || !finiteValue(frequencyTolerance)
     || referenceHeight <= VSMALL
     || angleTolerance < 0
     || calmFraction < 0
     || calmFraction > 1
     || frequencyTolerance < 0
    )
    {
        FatalIOErrorInFunction(climateDict)
            << "referenceHeight must be positive; directionTolerance and "
            << "frequencyTolerance must be non-negative; calmFraction "
            << "must be in [0,1]" << exit(FatalIOError);
    }

    forAll(directions, directionI)
    {
        if
        (
            mag(directions[directionI].referenceHeight - referenceHeight)
          > 1e-8*max(referenceHeight, scalar(1))
        )
        {
            FatalIOErrorInFunction(climateDict)
                << "Climate referenceHeight " << referenceHeight
                << " does not match referenceHeight "
                << directions[directionI].referenceHeight
                << " for windFrom " << directions[directionI].windFrom
                << ". Transform the climate data or CFD reference speed "
                << "to a common height first."
                << exit(FatalIOError);
        }
    }

    PtrList<dictionary> entries(climateDict.lookup("sectors"));
    List<ClimateSector> result(directions.size());
    boolList assigned(directions.size(), false);
    scalar totalFrequency = 0;

    forAll(entries, sectorI)
    {
        const dictionary& entry = entries[sectorI];
        ClimateSector sector;
        const scalar configuredWindFrom = entry.get<scalar>("windFrom");
        if (!finiteValue(configuredWindFrom))
        {
            FatalIOErrorInFunction(entry)
                << "windFrom must be finite"
                << exit(FatalIOError);
        }

        sector.windFrom = normalizeDirection(configuredWindFrom);
        sector.frequency = entry.get<scalar>("frequency");
        sector.scale = entry.get<scalar>("scale");
        sector.shape = entry.get<scalar>("shape");

        if
        (
           !finiteValue(sector.frequency)
         || !finiteValue(sector.scale)
         || !finiteValue(sector.shape)
         || sector.frequency < 0 || sector.frequency > 1
         || sector.scale <= VSMALL || sector.shape <= VSMALL
        )
        {
            FatalIOErrorInFunction(entry)
                << "frequency must be in [0,1], and scale/shape must be "
                << "positive for windFrom " << sector.windFrom
                << exit(FatalIOError);
        }

        label match = -1;
        scalar closest = GREAT;
        forAll(directions, directionI)
        {
            const scalar difference = directionDifference
            (
                sector.windFrom,
                directions[directionI].windFrom
            );
            if (difference < closest)
            {
                closest = difference;
                match = directionI;
            }
        }

        if (match < 0 || closest > angleTolerance)
        {
            FatalIOErrorInFunction(entry)
                << "No CFD direction matches climate windFrom "
                << sector.windFrom << " within " << angleTolerance
                << " degrees" << exit(FatalIOError);
        }

        if (assigned[match])
        {
            FatalIOErrorInFunction(entry)
                << "More than one climate sector matches CFD windFrom "
                << directions[match].windFrom << exit(FatalIOError);
        }

        result[match] = sector;
        assigned[match] = true;
        totalFrequency += sector.frequency;
    }

    forAll(assigned, directionI)
    {
        if (!assigned[directionI])
        {
            FatalIOErrorInFunction(climateDict)
                << "No climate sector supplied for CFD windFrom "
                << directions[directionI].windFrom
                << exit(FatalIOError);
        }
    }

    const scalar annualFrequency = totalFrequency + calmFraction;

    if (mag(annualFrequency - 1.0) > frequencyTolerance)
    {
        FatalIOErrorInFunction(climateDict)
            << "Climate sector frequencies sum to " << totalFrequency
            << " and calmFraction is " << calmFraction
            << ", giving " << annualFrequency << ". Their sum must equal "
            << "one within frequencyTolerance " << frequencyTolerance
            << exit(FatalIOError);
    }

    if (calmFraction > frequencyTolerance)
    {
        Info<< "Using explicit calm/non-exceeding annual fraction "
            << calmFraction << nl;
    }

    return result;
}


void writeResultField
(
    const word& fieldName,
    const scalarField& values,
    fvMesh& mesh
)
{
    if (values.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Cannot write field " << fieldName << ": received "
            << values.size() << " values for a mesh with " << mesh.nCells()
            << " cells" << exit(FatalError);
    }

    forAll(values, cellI)
    {
        if (!finiteValue(values[cellI]))
        {
            FatalErrorInFunction
                << "Cannot write field " << fieldName
                << ": non-finite value at cell " << cellI
                << exit(FatalError);
        }
    }

    volScalarField result
    {
        IOobject
        {
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        },
        mesh,
        dimensionedScalar(fieldName, dimless, 0),
        "zeroGradient"
    };

    result.primitiveFieldRef() = values;
    result.correctBoundaryConditions();

    Info<< "Writing " << fieldName << " to "
        << mesh.time().globalTimePath() << endl;

    if (!result.write())
    {
        FatalErrorInFunction
            << "Failed to write field " << fieldName
            << exit(FatalError);
    }
}


void validateMeshCoordinates
(
    const MeshSignature& reference,
    const fvMesh& mesh,
    const scalar tolerance,
    const scalar relativeVolumeTolerance,
    const fileName& caseDir
)
{
    if
    (
        mesh.nCells() != reference.cellCentres.size()
     || mesh.nPoints() != reference.nPoints
     || mesh.nFaces() != reference.nFaces
     || mesh.nInternalFaces() != reference.nInternalFaces
    )
    {
        FatalErrorInFunction
            << "Mesh sizes differ between the reference mesh and "
            << caseDir << nl
            << "Reference cells/points/faces/internalFaces: "
            << reference.cellCentres.size() << '/' << reference.nPoints
            << '/' << reference.nFaces << '/' << reference.nInternalFaces
            << nl << "Current cells/points/faces/internalFaces: "
            << mesh.nCells() << '/' << mesh.nPoints() << '/'
            << mesh.nFaces() << '/' << mesh.nInternalFaces() << nl
            << "Use matching meshes and coordinates."
            << exit(FatalError);
    }

    const vectorField& centres = mesh.cellCentres();
    const scalarField& volumes = mesh.V();
    scalar maxCentreDifference = 0;
    scalar maxRelativeVolumeDifference = 0;
    forAll(reference.cellCentres, cellI)
    {
        if
        (
           !finiteValue(centres[cellI])
         || !finiteValue(volumes[cellI])
         || volumes[cellI] <= VSMALL
        )
        {
            FatalErrorInFunction
                << "Mesh " << caseDir
                << " contains an invalid cell centre or volume at cell "
                << cellI << exit(FatalError);
        }

        maxCentreDifference = max
        (
            maxCentreDifference,
            mag(reference.cellCentres[cellI] - centres[cellI])
        );

        maxRelativeVolumeDifference = max
        (
            maxRelativeVolumeDifference,
            mag(reference.cellVolumes[cellI] - volumes[cellI])
           /max(reference.cellVolumes[cellI], volumes[cellI])
        );
    }

    if (maxCentreDifference > tolerance)
    {
        FatalErrorInFunction
            << "Cell coordinates/order differ for case " << caseDir
            << ". Maximum cell-centre difference is "
            << maxCentreDifference << " m, tolerance is " << tolerance
            << " m. Use matching meshes and coordinates."
            << exit(FatalError);
    }


    if (maxRelativeVolumeDifference > relativeVolumeTolerance)
    {
        FatalErrorInFunction
            << "Cell volumes/order differ for case " << caseDir
            << ". Maximum relative cell-volume difference is "
            << maxRelativeVolumeDifference << ", tolerance is "
            << relativeVolumeTolerance
            << ". Use matching meshes and coordinates."
            << exit(FatalError);
    }
}

} // End namespace Foam


int main(int argc, char *argv[])
{
    using namespace Foam;

    argList::addNote
    (
        "Evaluate pedestrian wind comfort from multiple directional cases"
    );
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary nenDict
    (
        IOobject
        (
            "NEN8100Dict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const fileName driverCase(runTime.globalPath());
    const word fieldName(nenDict.getOrDefault<word>("field", "U"));
    const word fieldPrefix
    (
        nenDict.getOrDefault<word>("fieldPrefix", "NEN8100")
    );

    const List<DirectionSpec> directions
    (
        readDirections(nenDict, driverCase)
    );

    const label expectedDirectionCount =
        nenDict.getOrDefault<label>("expectedDirectionCount", 12);
    const bool strictDirectionCount =
        nenDict.getOrDefault<bool>("strictDirectionCount", true);

    if (directions.size() != expectedDirectionCount)
    {
        if (strictDirectionCount)
        {
            FatalIOErrorInFunction(nenDict)
                << "Expected " << expectedDirectionCount
                << " directions but found " << directions.size()
                << exit(FatalIOError);
        }

        WarningInFunction
            << "Expected " << expectedDirectionCount
            << " directions for the configured assessment, but found "
            << directions.size() << nl;
    }

    const dictionary& surfaceDict = nenDict.subDict("evaluationSurface");
    const word surfaceType(surfaceDict.get<word>("type"));
    if (surfaceType != "horizontalPlane")
    {
        FatalIOErrorInFunction(surfaceDict)
            << "The initial implementation supports evaluationSurface type "
            << "horizontalPlane. Found '" << surfaceType << "'."
            << exit(FatalIOError);
    }

    const point groundPoint(surfaceDict.get<point>("groundPoint"));
    vector vertical(surfaceDict.get<vector>("verticalDirection"));
    const scalar pedestrianHeight =
        surfaceDict.getOrDefault<scalar>("pedestrianHeight", 1.75);

    if
    (
       !finiteValue(groundPoint)
     || !finiteValue(vertical)
     || !finiteValue(pedestrianHeight)
     || mag(vertical) <= VSMALL
     || pedestrianHeight < 0
    )
    {
        FatalIOErrorInFunction(surfaceDict)
            << "verticalDirection must be non-zero and pedestrianHeight "
            << "must be non-negative" << exit(FatalIOError);
    }

    vertical /= mag(vertical);
    const point slicePoint(groundPoint + pedestrianHeight*vertical);
    const bool triangulate =
        surfaceDict.getOrDefault<bool>("triangulate", true);
    dictionary sampledPlaneDict(surfaceDict);
    sampledPlaneDict.set("point", slicePoint);
    sampledPlaneDict.set("normal", vertical);
    sampledPlaneDict.set("triangulate", triangulate);
    const scalar geometryTolerance =
        nenDict.getOrDefault<scalar>("geometryTolerance", 1e-8);
    const scalar relativeVolumeTolerance =
        nenDict.getOrDefault<scalar>("relativeVolumeTolerance", 1e-10);

    if
    (
       !finiteValue(geometryTolerance)
     || !finiteValue(relativeVolumeTolerance)
     || geometryTolerance < 0
     || relativeVolumeTolerance < 0
    )
    {
        FatalIOErrorInFunction(nenDict)
            << "geometryTolerance and relativeVolumeTolerance must be "
            << "finite and non-negative"
            << exit(FatalIOError);
    }

    Info<< nl
        << "NEN 8100 pedestrian wind assessment" << nl
        << "  directions       : " << directions.size() << nl
        << "  pedestrian plane : point " << slicePoint
        << ", normal " << vertical << nl
        << "  output prefix    : " << fieldPrefix << nl
        << "  evaluation       : cell-centre velocity; the plane marks "
        << "intersected cells" << nl
        << endl;

    const bool hasClimate = nenDict.found("climate");
    List<ClimateSector> climate;
    if (hasClimate)
    {
        climate = readClimate(nenDict.subDict("climate"), directions);
    }

    const bool writeAmplification =
        nenDict.getOrDefault<bool>("writeAmplification", false);
    if (!hasClimate && !writeAmplification)
    {
        FatalIOErrorInFunction(nenDict)
            << "No climate dictionary is present and writeAmplification is "
            << "false, so there are no result fields to write."
            << exit(FatalIOError);
    }

    List<scalarField> amplification(directions.size());
    autoPtr<MeshSignature> referenceMesh;
    scalarField p5Fraction;
    scalarField p15Fraction;

    forAll(directions, directionI)
    {
        const DirectionSpec& spec = directions[directionI];

        Info<< "Sampling windFrom " << spec.windFrom << " deg from "
            << spec.caseDir << endl;

        Time caseTime
        (
            Time::controlDictName,
            spec.caseDir.path(),
            spec.caseDir.name(),
            false,
            true
        );
        selectCaseTime(caseTime, spec.timeSelection);

        fvMesh mesh
        (
            IOobject
            (
                spec.region,
                caseTime.timeName(),
                caseTime,
                IOobject::MUST_READ
            ),
            false
        );
        mesh.init(true);

        volVectorField U
        (
            IOobject
            (
                fieldName,
                caseTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        if (U.dimensions() != dimVelocity)
        {
            FatalErrorInFunction
                << "Field " << fieldName << " in " << spec.caseDir
                << " has dimensions " << U.dimensions()
                << ", expected velocity dimensions " << dimVelocity
                << exit(FatalError);
        }

        if (directionI == 0)
        {
            referenceMesh.reset(new MeshSignature(mesh));
            if (hasClimate)
            {
                p5Fraction.setSize(mesh.nCells(), 0);
                p15Fraction.setSize(mesh.nCells(), 0);
            }
        }
        else
        {
            validateMeshCoordinates
            (
                referenceMesh(),
                mesh,
                geometryTolerance,
                relativeVolumeTolerance,
                spec.caseDir
            );
        }

        scalarField gamma(mesh.nCells(), 0);
        const vectorField& velocity = U.primitiveField();
        forAll(velocity, cellI)
        {
            if (!finiteValue(velocity[cellI]))
            {
                FatalErrorInFunction
                    << "Field " << fieldName << " in " << spec.caseDir
                    << " contains a non-finite value at cell " << cellI
                    << exit(FatalError);
            }

            gamma[cellI] = mag(velocity[cellI])/spec.Uref;

            if (!finiteValue(gamma[cellI]))
            {
                FatalErrorInFunction
                    << "Non-finite velocity amplification in "
                    << spec.caseDir << " at cell " << cellI
                    << exit(FatalError);
            }

            if (hasClimate)
            {
                const ClimateSector& sector = climate[directionI];
                p5Fraction[cellI] += sector.frequency*weibullExceedance
                (
                    gamma[cellI],
                    5.0,
                    sector.scale,
                    sector.shape
                );
                p15Fraction[cellI] += sector.frequency*weibullExceedance
                (
                    gamma[cellI],
                    15.0,
                    sector.scale,
                    sector.shape
                );
            }
        }

        if (writeAmplification)
        {
            amplification[directionI].transfer(gamma);
        }
    }

    scalarField p5Percent(referenceMesh().cellCentres.size(), 0);
    scalarField p15Percent(referenceMesh().cellCentres.size(), 0);
    scalarField comfort(referenceMesh().cellCentres.size(), 0);
    scalarField traversing(referenceMesh().cellCentres.size(), 0);
    scalarField strolling(referenceMesh().cellCentres.size(), 0);
    scalarField sitting(referenceMesh().cellCentres.size(), 0);
    scalarField danger(referenceMesh().cellCentres.size(), 0);

    if (hasClimate)
    {
        forAll(referenceMesh().cellCentres, cellI)
        {
            p5Percent[cellI] = 100.0*p5Fraction[cellI];
            p15Percent[cellI] = 100.0*p15Fraction[cellI];
            const label cls = comfortClass(p5Percent[cellI]);
            comfort[cellI] = cls;
            traversing[cellI] = traversingRating(cls);
            strolling[cellI] = strollingRating(cls);
            sitting[cellI] = sittingRating(cls);
            danger[cellI] = dangerClass(p15Percent[cellI]);
        }
    }
    else
    {
        WarningInFunction
            << "No climate dictionary found. Writing directional "
            << "amplification factors only; no NEN classes are calculated."
            << nl;
    }

    fileName outputCase(directions[0].caseDir);
    if (nenDict.found("outputCase"))
    {
        outputCase = resolveCasePath
        (
            driverCase,
            nenDict.get<fileName>("outputCase")
        );
    }
    const word outputTimeSelection
    (
        nenDict.getOrDefault<word>("outputTime", "latest")
    );
    const word outputRegion
    (
        nenDict.getOrDefault<word>("outputRegion", directions[0].region)
    );

    if (!isDir(outputCase))
    {
        FatalIOErrorInFunction(nenDict)
            << "Output case directory does not exist: " << outputCase
            << exit(FatalIOError);
    }

    Time outputTime
    (
        Time::controlDictName,
        outputCase.path(),
        outputCase.name(),
        false,
        true
    );
    selectCaseTime(outputTime, outputTimeSelection);

    fvMesh outputMesh
    (
        IOobject
        (
            outputRegion,
            outputTime.timeName(),
            outputTime,
            IOobject::MUST_READ
        ),
        false
    );
    outputMesh.init(true);

    validateMeshCoordinates
    (
        referenceMesh(),
        outputMesh,
        geometryTolerance,
        relativeVolumeTolerance,
        outputCase
    );

    sampledPlane pedestrianPlane
    (
        "NEN8100PedestrianPlane",
        outputMesh,
        sampledPlaneDict
    );
    pedestrianPlane.update();

    if (pedestrianPlane.faces().empty())
    {
        FatalErrorInFunction
            << "The pedestrian plane does not intersect output case "
            << outputCase << exit(FatalError);
    }

    scalarField pedestrianMask(outputMesh.nCells(), 0);
    const labelList& cutCells = pedestrianPlane.meshCells();
    scalar maxCellCentreOffset = 0;
    forAll(cutCells, cutI)
    {
        if (cutCells[cutI] >= 0)
        {
            pedestrianMask[cutCells[cutI]] = 1;
            maxCellCentreOffset = max
            (
                maxCellCentreOffset,
                mag
                (
                    (outputMesh.cellCentres()[cutCells[cutI]] - slicePoint)
                  & vertical
                )
            );
        }
    }

    Info<< "Maximum normal distance from the pedestrian plane to a marked "
        << "cell centre: " << maxCellCentreOffset << " m" << nl;

    writeResultField
    (
        fieldPrefix + "PedestrianMask",
        pedestrianMask,
        outputMesh
    );

    if (writeAmplification)
    {
        forAll(directions, directionI)
        {
            writeResultField
            (
                directionFieldName(fieldPrefix, directions[directionI].windFrom),
                amplification[directionI],
                outputMesh
            );
        }
    }

    if (hasClimate)
    {
        writeResultField(fieldPrefix + "P5Percent", p5Percent, outputMesh);
        writeResultField(fieldPrefix + "ComfortClass", comfort, outputMesh);
        writeResultField
        (
            fieldPrefix + "TraversingRating",
            traversing,
            outputMesh
        );
        writeResultField
        (
            fieldPrefix + "StrollingRating",
            strolling,
            outputMesh
        );
        writeResultField
        (
            fieldPrefix + "SittingRating",
            sitting,
            outputMesh
        );
        writeResultField(fieldPrefix + "P15Percent", p15Percent, outputMesh);
        writeResultField(fieldPrefix + "DangerClass", danger, outputMesh);
    }

    Info<< nl
        << "Wrote OpenFOAM result fields to " << outputCase
        << '/' << outputTime.timeName() << nl
        << "Use " << fieldPrefix << "PedestrianMask = 1 to select the "
        << "configured pedestrian-height cells in ParaView." << nl;
    if (hasClimate)
    {
        Info<< "Class encoding: " << fieldPrefix
            << "ComfortClass 1=A ... 5=E" << nl
            << "Activity encoding: 0=good, 1=moderate, 2=poor" << nl
            << "Danger encoding: 0=no danger, 1=limited risk, 2=dangerous"
            << nl;
    }
    Info<< endl;

    return 0;
}

// ************************************************************************* //
