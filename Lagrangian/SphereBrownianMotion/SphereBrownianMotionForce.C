/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SphereBrownianMotionForce.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "demandDrivenData.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SphereBrownianMotionForce<CloudType>::erfInv(const scalar y) const
{
    const scalar a = 0.147;
    scalar k = 2.0/(mathematical::pi*a) +  0.5*log(1.0 - y*y);
    scalar h = log(1.0 - y*y)/a;
    scalar x = sqrt(-k + sqrt(k*k - h));

    if (y < 0.0)
    {
        return -x;
    }
    else
    {
        return x;
    }
}

template<class CloudType>
Foam::SphereBrownianMotionForce<CloudType>::SphereBrownianMotionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    rndGen_(owner.rndGen())
{}


template<class CloudType>
Foam::SphereBrownianMotionForce<CloudType>::SphereBrownianMotionForce
(
    const SphereBrownianMotionForce& sbmf
)
:
    ParticleForce<CloudType>(sbmf),
    rndGen_(sbmf.rndGen_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SphereBrownianMotionForce<CloudType>::~SphereBrownianMotionForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::SphereBrownianMotionForce<CloudType>::calcCoupled
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
    const scalar dp = p.d();
    const scalar Tc = td.Tc();
    scalar f = 0;
    f=sqrt(6*mathematical::pi*muc*kb*Tc*dp/dt);

    // To generate a cubic distribution (3 independent directions) :
    const scalar sqrt2 = sqrt(2.0);
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        const scalar x = rndGen_.sample01<scalar>();
        const scalar eta = sqrt2*erfInv(2*x - 1.0);
        value.Su()[dir] = f*eta;
    }

    return value;
}


// ************************************************************************* //
