/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "fixedHumidityFvPatchScalarField.H"

class heHumidityRhoThermo;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    mode_("relative"),
    method_("buck"),
    value_(0.0),
    methodName_
    (
        IOobject
        (
            "methodName",
            p.boundaryMesh().mesh().time().timeName(),
            p.boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1
    )
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    mode_(ptf.mode_),
    method_(ptf.method_),
    value_(ptf.value_),
    methodName_(ptf.methodName_)
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    mode_(dict.lookupOrDefault<word>("mode", "relative")),
    method_(dict.lookupOrDefault<word>("method", "buck")),
    value_(readScalar(dict.lookup("humidity"))),
    methodName_
    (
        IOobject
        (
            "methodName",
            p.boundaryMesh().mesh().time().timeName(),
            p.boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1
    )
{
    if (mode_ == "absolute")
    {
        NotImplemented;
    }
    else if (mode_ != "relative")
    {
        FatalErrorInFunction
            << "The specified type is not supported '"
            << mode_ << "'. Supported are 'relative' or 'absolute'"
            << exit(FatalError);
    }

    if (method_ != "buck" && method_ != "magnus")
    {
        FatalErrorInFunction
            << "The specified method to calculate the saturation pressure is "
            << "not supported: " << method_ << ". Supported methods are "
            << "'buck' and 'magnus'."
            << exit(FatalError);
    }

    methodName_[0] = method_;
}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    mode_(tppsf.mode_),
    method_(tppsf.method_),
    value_(tppsf.value_),
    methodName_(tppsf.methodName_)
{}


Foam::fixedHumidityFvPatchScalarField::
fixedHumidityFvPatchScalarField
(
    const fixedHumidityFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    mode_(tppsf.mode_),
    method_(tppsf.method_),
    value_(tppsf.value_),
    methodName_(tppsf.methodName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedHumidityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const label patchi = patch().index();

    const scalarField specificHumidity = calcSpecificHumidity(thermo, patchi);

    //const scalarField& pw = thermo.p().boundaryField()[patchi];
    //fvPatchScalarField& Tw =
    //    const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    //Tw.evaluate();
    operator==(specificHumidity);

    fixedValueFvPatchScalarField::updateCoeffs();
}


const Foam::scalarField Foam::fixedHumidityFvPatchScalarField::calcSpecificHumidity
(
    const basicThermo& thermo,
    const label patchi
)
{
    //- a) Calc saturation pressure
    const scalarField& pfT = thermo.T().boundaryField()[patchi];
    const scalarField& pfp = thermo.p().boundaryField()[patchi];
    const scalarField theta = pfT - 273.15; //(Tpatch.size(), scalar(0));

    scalarField pSatH2O(pfT.size(), scalar(0));

    //  Standard method not as accurate as buck
    if (method_ == "magnus")
    {
        scalar pre1 = 611.2;
        scalar pre2 = 17.62;
        scalar value1 = 243.12;

        forAll(pSatH2O, facei)
        {
            pSatH2O[facei] =
                pre1*exp((pre2*(theta[facei]))/(value1+theta[facei]));
        }
    }
    //  Buck formula [1996]
    //  Valid between 0 to 100 degC and 1013.25 hPa
    //  Very accurate between 0 degC and 50 degC
    else if (method_ == "buck")
    {
        scalar pre1 = 611.21;
        scalar value1 = 18.678;
        scalar value2 = 234.5;
        scalar value3 = 257.14;

        forAll(pSatH2O, facei)
        {
            scalar TdC = theta[facei];

            pSatH2O[facei] =
                pre1*exp(((value1-TdC/value2)*TdC/(value3+TdC)));
        }
    }

    //- b) Calc partial pressure of water
    scalarField partialPressureH2O(pfT.size(), scalar(0));
    partialPressureH2O = value_ * pSatH2O;

    //- c) Calc density of water [kg/m^3]
    scalarField rhoWater(pfT.size(), scalar(0));
    {
        scalar RH2O = 461.51;
        rhoWater = partialPressureH2O/(RH2O*pfT);
    }

    //- d) Calc density of dry air [kg/m^3]
    scalarField rhoDryAir(pfT.size(), scalar(0));
    {
        scalar RdryAir = 287.058;
        rhoDryAir = (pfp - partialPressureH2O)/(RdryAir*pfT);
    }

    //- e) Calculate the specific humidity [kg/kg]
    scalarField specificHumidity(pfT.size(), scalar(0));

    return rhoWater/(rhoWater+rhoDryAir);
}


void Foam::fixedHumidityFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("mode") << mode_ << token::END_STATEMENT << nl;
    os.writeKeyword("method") << method_ << token::END_STATEMENT << nl;
    os.writeKeyword("humidity") << value_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedHumidityFvPatchScalarField
    );
}

// ************************************************************************* //
