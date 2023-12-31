/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HeatTransferModel<CloudType>::HeatTransferModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    BirdCorrection_(false)
{}


template<class CloudType>
Foam::HeatTransferModel<CloudType>::HeatTransferModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    BirdCorrection_(coeffDict_.lookup("BirdCorrection"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HeatTransferModel<CloudType>::~HeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::HeatTransferModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::HeatTransferModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::HeatTransferModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
const Foam::Switch& Foam::HeatTransferModel<CloudType>::BirdCorrection() const
{
    return BirdCorrection_;
}


template<class CloudType>
Foam::scalar Foam::HeatTransferModel<CloudType>::htc
(
    const scalar dp,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW
) const
{
    const scalar Nu = this->Nu(Re, Pr);

    scalar htc = Nu*kappa/dp;

    if (BirdCorrection_ && (mag(htc) > ROOTVSMALL) && (mag(NCpW) > ROOTVSMALL))
    {
        const scalar phit = min(NCpW/htc, 50);
        if (phit > 0.001)
        {
            htc *= phit/(exp(phit) - 1.0);
        }
    }

    return htc;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewHeatTransferModel.C"


// ************************************************************************* //

