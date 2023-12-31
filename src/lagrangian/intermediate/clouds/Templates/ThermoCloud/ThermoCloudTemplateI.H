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

#include "radiationConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::ThermoCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::basicThermo&
Foam::ThermoCloud<ParcelType>::carrierThermo() const
{
    return carrierThermo_;
}


template<class ParcelType>
inline Foam::basicThermo&
Foam::ThermoCloud<ParcelType>::carrierThermo()
{
    return carrierThermo_;
}


template<class ParcelType>
inline const Foam::HeatTransferModel<Foam::ThermoCloud<ParcelType> >&
Foam::ThermoCloud<ParcelType>::heatTransfer() const
{
    return heatTransferModel_;
}


template<class ParcelType>
inline const Foam::scalarIntegrationScheme&
Foam::ThermoCloud<ParcelType>::TIntegrator() const
{
    return TIntegrator_;
}


template<class ParcelType>
inline bool Foam::ThermoCloud<ParcelType>::radiation() const
{
    return radiation_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<ParcelType>::hsTrans()
{
    return hsTrans_;
}


template<class ParcelType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ThermoCloud<ParcelType>::Sh() const
{
    tmp<DimensionedField<scalar, volMesh> > tSh
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "Sh",
                this->db().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            hsTrans_/(this->mesh().V()*this->db().time().deltaT())
        )
    );

    return tSh;
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<ParcelType>::Ep() const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "radiation::Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    // Need to check if coupled as field is created on-the-fly
    if (radiation_ && this->coupled())
    {
        scalarField& Ep = tEp().internalField();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();

        forAllConstIter(typename ThermoCloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
            const label cellI = p.cell();
            Ep[cellI] += p.nParticle()*p.areaP()*pow4(p.T());
        }

        Ep *= epsilon*radiation::sigmaSB.value()/V;
    }

    return tEp;
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<ParcelType>::ap() const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "radiation::ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    // Need to check if coupled as field is created on-the-fly
    if (radiation_ && this->coupled())
    {
        scalarField& ap = tap().internalField();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();

        forAllConstIter(typename ThermoCloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
            const label cellI = p.cell();
            ap[cellI] += p.nParticle()*p.areaP();
        }

        ap *= epsilon/V;
    }

    return tap;
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<ParcelType>::sigmap() const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "radiation::sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    // Need to check if coupled as field is created on-the-fly
    if (radiation_ && this->coupled())
    {
        scalarField& sigmap = tsigmap().internalField();

        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalar f = constProps_.f0();

        forAllConstIter(typename ThermoCloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
            const label cellI = p.cell();
            sigmap[cellI] += p.nParticle()*p.areaP();
        }

        sigmap *= (1.0 - f)*(1.0 - epsilon)/V;
    }

    return tsigmap;
}


// ************************************************************************* //
