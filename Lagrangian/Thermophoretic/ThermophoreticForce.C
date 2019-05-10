/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "ThermophoreticForce.H"
#include "fundamentalConstants.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermophoreticForce<CloudType>::ThermophoreticForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    TName_(this->coeffs().template lookupOrDefault<word>("T", "T")),
    gradTInterpPtr_(nullptr),
    ST_(readScalar(this->coeffs().lookup("ST")))
{}


template<class CloudType>
Foam::ThermophoreticForce<CloudType>::ThermophoreticForce
(
    const ThermophoreticForce& tpf
)
:
    ParticleForce<CloudType>(tpf),
    TName_(tpf.TName_),
    gradTInterpPtr_(nullptr),
    ST_(tpf.ST_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermophoreticForce<CloudType>::~ThermophoreticForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermophoreticForce<CloudType>::cacheFields(const bool store)
{
    static word fName("gradT");

    bool fieldExists = this->mesh().template foundObject<volVectorField>(fName);

    if (store)
    {
        if (!fieldExists)
        {
            const volScalarField& Tc = this->mesh().template
                lookupObject<volScalarField>(TName_);

            volVectorField* gradTPtr = new volVectorField
            (
                fName,
                (fvc::grad(Tc))
            );

            gradTPtr->store();
        }

        const volVectorField& gradT = this->mesh().template
            lookupObject<volVectorField>(fName);

        gradTInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                gradT
            ).ptr()
        );
    }
    else
    {
        gradTInterpPtr_.clear();

        if (fieldExists)
        {
            volVectorField& gradT = 
                this->mesh().template lookupObjectRef<volVectorField>(fName);

            gradT.checkOut();
        }
    }

}


template<class CloudType>
Foam::forceSuSp Foam::ThermophoreticForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    const scalar kb = physicoChemical::k.value();
    vector gradT = 
        gradTInterp().interpolate(p.coordinates(), p.currentTetIndices());
    //scalar D=kb*td.Tc/(3*mathematical::pi*muc*dp);
    //value.Su() = -3.0*mathematical::pi*muc*dp*D*ST()*gradT;
    value.Su() = -1.0*kb*td.Tc()*ST_*gradT;

    return value;    
}

template<class CloudType>
Foam::scalar Foam::ThermophoreticForce<CloudType>::massAdd
(
    const typename CloudType::parcelType&,
    const typename CloudType::parcelType::trackingData& td,
    const scalar
) const
{
    return 0.0;
}

// ************************************************************************* //
