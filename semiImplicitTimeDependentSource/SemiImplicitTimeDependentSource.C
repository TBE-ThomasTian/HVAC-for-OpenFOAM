/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "SemiImplicitTimeDependentSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeModeTypeNames_
{
    "absolute", "specific"
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
typename Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeModeType
Foam::fv::SemiImplicitTimeDependentSource<Type>::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


template<class Type>
Foam::word Foam::fv::SemiImplicitTimeDependentSource<Type>::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}


template<class Type>
void Foam::fv::SemiImplicitTimeDependentSource<Type>::setFieldData(const dictionary& dict)
{
    label count = dict.size();

    fieldNames_.resize(count);
    injectionRate_.resize(count);
    applied_.resize(count, false);

    count = 0;
    for (const entry& dEntry : dict)
    {
        fieldNames_[count] = dEntry.keyword();
        dEntry.readEntry(injectionRate_[count]);

        ++count;
    }

    // Set volume normalisation
    if (volumeMode_ == vmAbsolute)
    {
        VDash_ = V_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SemiImplicitTimeDependentSource<Type>::SemiImplicitTimeDependentSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
    injectionRate_(),
	timeDependence_("constant/heatReleaseRate")//JN: read in HRR file with the HRR
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitTimeDependentSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitTimeDependentSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    typename GeometricField<Type, fvPatchField, volMesh>::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>(eqn.dimensions()/dimVolume, Zero),
        false
    );

	//JN: here we multiply the explicit source S_u from fvOptions with HRR
	//JN: always use 1 for S_u and 0 for S_p if you define a time dependent HRR
    UIndirectList<Type>(Su, cells_) = timeDependence_(mesh_.time().value())*injectionRate_[fieldi].first()/VDash_;

    volScalarField::Internal Sp
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>(Su.dimensions()/psi.dimensions(), Zero),
        false
    );

    UIndirectList<scalar>(Sp, cells_) = injectionRate_[fieldi].second()/VDash_;

    eqn += Su + fvm::SuSp(Sp, psi);
}


template<class Type>
void Foam::fv::SemiImplicitTimeDependentSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitTimeDependentSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldi);
}



template<class Type>
bool Foam::fv::SemiImplicitTimeDependentSource<Type>::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        volumeMode_ = wordToVolumeModeType(coeffs_.get<word>("volumeMode"));
        setFieldData(coeffs_.subDict("injectionRateSuSp"));

        return true;
    }

    return false;
}


// ************************************************************************* //
