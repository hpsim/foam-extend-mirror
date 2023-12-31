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

#include "particleForces.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleForces::particleForces
(
    const fvMesh& mesh,
    const dictionary& dict,
    const vector& g
)
:
    mesh_(mesh),
    dict_(dict.subDict("particleForces")),
    g_(g),
    gradUPtr_(nullptr),
    gravity_(dict_.lookup("gravity")),
    virtualMass_(dict_.lookup("virtualMass")),
    Cvm_(0.0),
    pressureGradient_(dict_.lookup("pressureGradient")),
    UName_(dict_.lookupOrDefault<word>("U", "U"))
{
    if (virtualMass_)
    {
        dict_.lookup("Cvm") >> Cvm_;
    }
}


Foam::particleForces::particleForces(const particleForces& f)
:
    mesh_(f.mesh_),
    dict_(f.dict_),
    g_(f.g_),
    gradUPtr_(f.gradUPtr_),
    gravity_(f.gravity_),
    virtualMass_(f.virtualMass_),
    Cvm_(f.Cvm_),
    pressureGradient_(f.pressureGradient_),
    UName_(f.UName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleForces::~particleForces()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::particleForces::dict() const
{
    return dict_;
}


const Foam::vector& Foam::particleForces::g() const
{
    return g_;
}


Foam::Switch Foam::particleForces::gravity() const
{
    return gravity_;
}


Foam::Switch Foam::particleForces::virtualMass() const
{
    return virtualMass_;
}


Foam::Switch Foam::particleForces::pressureGradient() const
{
    return pressureGradient_;
}


const Foam::word& Foam::particleForces::UName() const
{
    return UName_;
}


void Foam::particleForces::cacheFields(const bool store)
{
    if (store && pressureGradient_)
    {
        const volVectorField U = mesh_.lookupObject<volVectorField>(UName_);
        gradUPtr_ = fvc::grad(U).ptr();
    }
    else
    {
        if (gradUPtr_)
        {
            delete gradUPtr_;
            gradUPtr_ = nullptr;
        }
    }
}


Foam::vector Foam::particleForces::calcCoupled
(
    const label cellI,
    const scalar dt,
    const scalar rhoc,
    const scalar rho,
    const vector& Uc,
    const vector& U
) const
{
    vector Ftot = vector::zero;

    // Virtual mass force
    if (virtualMass_)
    {
        notImplemented
        (
            "Foam::particleForces::calcCoupled(...) - virtual mass force"
        );
//        Ftot += Cvm_*rhoc/rho*d(Uc - U)/dt;
    }

    // Pressure gradient force
    if (pressureGradient_)
    {
        const volTensorField& gradU = *gradUPtr_;
        Ftot += rhoc/rho*(U & gradU[cellI]);
    }

    return Ftot;
}


Foam::vector Foam::particleForces::calcNonCoupled
(
    const label cellI,
    const scalar dt,
    const scalar rhoc,
    const scalar rho,
    const vector& Uc,
    const vector& U
) const
{
    vector Ftot = vector::zero;

    // Gravity force
    if (gravity_)
    {
        Ftot += g_*(1.0 - rhoc/rho);
    }

    return Ftot;
}


// ************************************************************************* //

