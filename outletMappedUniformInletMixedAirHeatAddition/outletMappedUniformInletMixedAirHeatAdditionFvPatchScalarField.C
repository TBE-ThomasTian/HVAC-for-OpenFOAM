/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#include "outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

static Foam::word firstPresentKeyword
(
    const Foam::dictionary& dict,
    const Foam::word& primary,
    const Foam::word& alias
)
{
    if (dict.found(primary))
    {
        return primary;
    }

    if (dict.found(alias))
    {
        return alias;
    }

    return primary;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Qptr_(nullptr),
    additionalMassFlowRatePtr_(nullptr),
    additionalTemperaturePtr_(nullptr),
    outletPatchName_(),
    phiName_("phi"),
    additionalCp_(-1),
    TMin_(0),
    TMax_(5000),
    debug_(false),
    debugTimeIndex_(-1)
{}


Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
(
    const outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Qptr_(ptf.Qptr_.clone()),
    additionalMassFlowRatePtr_(ptf.additionalMassFlowRatePtr_.clone()),
    additionalTemperaturePtr_(ptf.additionalTemperaturePtr_.clone()),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    additionalCp_(ptf.additionalCp_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_),
    debug_(ptf.debug_),
    debugTimeIndex_(ptf.debugTimeIndex_)
{}


Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Qptr_(Function1<scalar>::NewIfPresent("Q", dict, word::null, &db())),
    additionalMassFlowRatePtr_
    (
        Function1<scalar>::New
        (
            firstPresentKeyword
            (
                dict,
                "additionalMassFlowRate",
                "freshMassFlowRate"
            ),
            dict,
            &db()
        )
    ),
    additionalTemperaturePtr_
    (
        Function1<scalar>::New
        (
            firstPresentKeyword
            (
                dict,
                "additionalTemperature",
                "freshTemperature"
            ),
            dict,
            &db()
        )
    ),
    outletPatchName_(dict.get<word>("outletPatch")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    additionalCp_(dict.getOrDefault<scalar>("additionalCp", -1)),
    TMin_(dict.getOrDefault<scalar>("TMin", 0)),
    TMax_(dict.getOrDefault<scalar>("TMax", 5000)),
    debug_(dict.getOrDefault<Switch>("debug", false)),
    debugTimeIndex_(-1)
{}


Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
(
    const outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    Qptr_(ptf.Qptr_.clone()),
    additionalMassFlowRatePtr_(ptf.additionalMassFlowRatePtr_.clone()),
    additionalTemperaturePtr_(ptf.additionalTemperaturePtr_.clone()),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    additionalCp_(ptf.additionalCp_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_),
    debug_(ptf.debug_),
    debugTimeIndex_(ptf.debugTimeIndex_)
{}


Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
(
    const outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    Qptr_(ptf.Qptr_.clone()),
    additionalMassFlowRatePtr_(ptf.additionalMassFlowRatePtr_.clone()),
    additionalTemperaturePtr_(ptf.additionalTemperaturePtr_.clone()),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    additionalCp_(ptf.additionalCp_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_),
    debug_(ptf.debug_),
    debugTimeIndex_(ptf.debugTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    const volScalarField& T =
    (
        dynamic_cast<const volScalarField&>(this->internalField())
    );

    const fvPatch& fvp = patch();

    const label outletPatchID =
        fvp.patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find outlet patch " << outletPatchName_
            << abort(FatalError);
    }

    const fvPatch& outletPatch = fvp.boundaryMesh()[outletPatchID];

    const fvPatchField<scalar>& outletPatchField =
        T.boundaryField()[outletPatchID];

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField& outletPatchPhi = phi.boundaryField()[outletPatchID];

    const basicThermo& thermo =
        db().lookupObject<basicThermo>(basicThermo::dictName);

    const scalarField& pp = thermo.p().boundaryField()[outletPatchID];
    const scalarField& pT = thermo.T().boundaryField()[outletPatchID];
    const tmp<scalarField> tCpf(thermo.Cp(pp, pT, outletPatchID));
    const scalarField& Cpf = tCpf();

    scalarField massWeights(outletPatchPhi);

    if (phi.dimensions() == dimVolume/dimTime)
    {
        const tmp<scalarField> trhop(thermo.rho(outletPatchID));
        massWeights *= trhop();
    }
    else if (phi.dimensions() != dimMass/dimTime)
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    const scalar recircMassFlow = gSum(massWeights);

    scalar averageOutletTemperature =
        gWeightedAverage(outletPatch.magSf(), outletPatchField);

    scalar recircCp = gAverage(Cpf);

    if (recircMassFlow > SMALL)
    {
        averageOutletTemperature =
            gWeightedSum(massWeights, outletPatchField)/recircMassFlow;

        recircCp = gWeightedSum(massWeights, Cpf)/recircMassFlow;
    }

    const scalar Q = Qptr_ ? Qptr_->value(t) : scalar(0);

    const scalar additionalMassFlowRate =
        additionalMassFlowRatePtr_->value(t);

    if (additionalMassFlowRate < -SMALL)
    {
        FatalErrorInFunction
            << "additionalMassFlowRate must be non-negative. Got "
            << additionalMassFlowRate
            << nl << exit(FatalError);
    }

    const scalar additionalTemperature =
        additionalTemperaturePtr_->value(t);

    scalar additionalCp = additionalCp_;

    if (additionalCp <= SMALL)
    {
        const scalarField& inletP =
            thermo.p().boundaryField()[this->patch().index()];

        const scalarField inletAdditionalT(this->patch().size(), additionalTemperature);

        const tmp<scalarField> tCpAdditional
        (
            thermo.Cp(inletP, inletAdditionalT, this->patch().index())
        );

        additionalCp = gAverage(tCpAdditional());
    }

    const scalar recircHeatCapacityRate = recircMassFlow*recircCp;
    const scalar additionalHeatCapacityRate =
        max(additionalMassFlowRate, scalar(0))*additionalCp;

    const scalar totalHeatCapacityRate =
        recircHeatCapacityRate + additionalHeatCapacityRate;

    scalar inletTemperature = averageOutletTemperature;

    if (totalHeatCapacityRate > SMALL)
    {
        inletTemperature =
        (
            recircHeatCapacityRate*averageOutletTemperature
          + additionalHeatCapacityRate*additionalTemperature
          + Q
        )/totalHeatCapacityRate;
    }

    operator==(clamp(inletTemperature, TMin_, TMax_));

    if (debug_)
    {
        const label ti = db().time().timeIndex();
        if (debugTimeIndex_ != ti)
        {
            debugTimeIndex_ = ti;

            Info<< typeName << " " << patch().name()
                << " T_out=" << averageOutletTemperature
                << " mdot_recirc=" << recircMassFlow
                << " Cp_recirc=" << recircCp
                << " mdot_add=" << additionalMassFlowRate
                << " T_add=" << additionalTemperature
                << " Cp_add=" << additionalCp
                << " Q=" << Q
                << " T_in=" << clamp(inletTemperature, TMin_, TMax_)
                << endl;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("outletPatch", outletPatchName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);

    if (Qptr_)
    {
        Qptr_->writeData(os);
    }
    else
    {
        os.writeEntry("Q", scalar(0));
    }

    additionalMassFlowRatePtr_->writeData(os);
    additionalTemperaturePtr_->writeData(os);

    if (additionalCp_ > SMALL)
    {
        os.writeEntry("additionalCp", additionalCp_);
    }

    os.writeEntry("TMin", TMin_);
    os.writeEntry("TMax", TMax_);

    if (debug_)
    {
        os.writeEntry("debug", debug_);
    }

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletMappedUniformInletMixedAirHeatAdditionFvPatchScalarField
    );
}


// ************************************************************************* //
