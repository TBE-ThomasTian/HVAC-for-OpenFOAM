/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "heHumidityRhoThermo.H"
#include "fvMatricesFwd.H"
#include "fvCFD.H"
#include "bound.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::heHumidityRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::heHumidityRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName, dictName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::~heHumidityRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    Info<< "   Solve transport equation for specific humidity\n";
    specificHumidityTransport();

    Info<< "   Calculate the saturation pressure of water\n";
    pSatH2O();

    Info<< "   Calculate the partial pressure of water\n";
    partialPressureH2O();

    Info<< "   Calculate the relative humidity\n";
    relHumidity();

    Info<< "   Calculate the water content\n";
    waterContent();

    Info<< "   Accounting for density change based on humidity\n";
    densityChange();

    //- Keep physical bounds
    limit();

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::pSatH2O()
{
    const volScalarField theta =
        (this->T_
      - dimensionedScalar("Kelvin", dimensionSet(0,0,0,1,0,0,0), scalar(273.15)))
      * dimensionedScalar("perK", dimensionSet(0,0,0,-1,0,0,0), scalar(1));

    //- Magnus formulation
    //  Valid between -50 to 100 degC and 1013.25 hPa
    if (this->method_ == "magnus")
    {
        dimensionedScalar pre1
        (
            "pre1",
            dimensionSet(1,-1,-2,0,0,0,0),
            scalar(611.2)
        );

        dimensionedScalar pre2
        (
            "pre2",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(17.62)
        );

        dimensionedScalar value1
        (
            "value1",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(243.12)
        );

        this->pSatH2O_ = pre1*exp((pre2*(theta))/(value1+theta));
    }
    //  Buck formula [1996]
    //  Valid between 0 to 100 degC and 1013.25 hPa
    //  Very accurate between 0 degC and 50 degC
    else if (this->method_ == "buck")
    {
        dimensionedScalar pre1
        (
            "pre1",
            dimensionSet(1,-1,-2,0,0,0,0),
            scalar(611.21)
        );

        dimensionedScalar value1
        (
            "value1",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(18.678)
        );

        dimensionedScalar value2
        (
            "value2",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(234.5)
        );

        dimensionedScalar value3
        (
            "value3",
            dimensionSet(0,0,0,0,0,0,0),
            scalar(257.14)
        );

        this->pSatH2O_ =
            pre1*exp(((value1-(theta/value2))*theta/(value3+theta)));
    }
    else
    {
        FatalErrorInFunction
            << "The specified method to calculate the saturation pressure is "
            << "not supported: " << this->method_ << ". Supported methods are "
            << "'buck' and 'magnus'."
            << exit(FatalError);
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::partialPressureH2O()
{
    const dimensionedScalar RSpecificH2O
    (
        "gasConstantH2O",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    const dimensionedScalar RSpecificDryAir
    (
        "gasConstantDryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    const volScalarField& p = this->p_;
    const volScalarField& T = this->T_;
    const volScalarField& sH = this->specificHumidity_;

    volScalarField& pPH2O = this->partialPressureH2O_;

    pPH2O =
        (-p/RSpecificDryAir/T)
      * pow
        (
            (-1./(RSpecificDryAir*T))
          + ((1-pow(sH+VSMALL, -1))* (1/(RSpecificH2O*T))),
           -1
        );
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::relHumidity()
{
    volScalarField& relHum = this->relHum_;

    relHum = this->partialPressureH2O_ / this->pSatH2O_;

    forAll(relHum, cellI)
    {
        if (relHum[cellI] > 1.)
        {
            WarningInFunction
                << "Humidity exeeds 100 percent, condensation occur in cell "
                << cellI << " (not implemented)" << endl;
            break;
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::waterContent()
{
    const volScalarField& pPH2O = this->partialPressureH2O_;
    const volScalarField& pSatH2O = this->pSatH2O_;
    const volScalarField& p = this->p_;
    const volScalarField& T= this->T_;

    const dimensionedScalar RSpecificH2O
    (
        "gasConstantH2O",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(461.51)
    );

    volScalarField& waterContent = this->waterContent_;
    waterContent = pPH2O / (RSpecificH2O * T);

    volScalarField& maxWaterContent = this->maxWaterContent_;
    maxWaterContent = pSatH2O / (RSpecificH2O * T);

    const dimensionedScalar RSpecificDryAir
    (
        "gasConstantDryAir",
        dimensionSet(0,2,-2,-1,0,0,0),
        scalar(287.058)
    );

    volScalarField& maxSpecificHumidity = this->maxSpecificHumidity_;
    maxSpecificHumidity =
        maxWaterContent / ((p-pSatH2O)/(RSpecificDryAir*T) + maxWaterContent);
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::
specificHumidityTransport()
{
    volScalarField& specHum = this->specificHumidity_;
    volScalarField& muEff = this->muEff_;

    const volScalarField& mu = this->mu_;

    //- Old density
    const volScalarField& rho =
        this->db().objectRegistry::lookupObject<volScalarField>("rho");

    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>("phi");

    const IOdictionary& turbProp =
        this->db().objectRegistry
        ::lookupObject<IOdictionary>("turbulenceProperties");

    const word turbulenceMode = word(turbProp.lookup("simulationType"));

    if (turbulenceMode == "RAS")
    {
        const volScalarField& nut =
            this->db().objectRegistry::lookupObject<volScalarField>("nut");

        muEff = rho*nut + mu;
    }
    else
    {
        muEff = mu;
    }

    fvScalarMatrix specHumEqn
    (
        fvm::ddt(rho, specHum)
      + fvm::div(phi, specHum)
     ==
        fvm::laplacian(muEff, specHum)
    );

    specHumEqn.relax();
    specHumEqn.solve();

    //- To keep physical range
    //  Defined between 0 and max water content based on saturation pressure
    //  Maximum is not limited yet
    bound
    (
        specHum,
        dimensionedScalar("tmp", dimensionSet(0,0,0,0,0,0,0), scalar(0))
    );
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::densityChange()
{
    this->rho_ += this->waterContent_;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heHumidityRhoThermo<BasicPsiThermo, MixtureType>::limit()
{
    volScalarField& specHum = this->specificHumidity_;
    const volScalarField& maxSpecHum = this->maxSpecificHumidity_;

    label min = 0;
    label max = 0;

    forAll(specHum, cellI)
    {
        if (specHum[cellI] < 0)
        {
            specHum[cellI] = 0;
            min++;
        }
        else if (specHum[cellI] > maxSpecHum[cellI])
        {
            specHum[cellI] = maxSpecHum[cellI];
            max++;
        }
    }

    if (min > 0)
    {
        Info<< "Corrected " << min << " cells which were lower than 0\n";
    }
    if (max > 0)
    {
        Info<< "Corrected " << max << " cells which were higher than max\n";
    }
}

// ************************************************************************* //
