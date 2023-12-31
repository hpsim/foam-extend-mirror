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

#include "JohnsonJacksonFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(JohnsonJacksonFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        JohnsonJacksonFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonFrictionalStress::JohnsonJacksonFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::JohnsonJacksonFrictionalStress::~JohnsonJacksonFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::JohnsonJacksonFrictionalStress::
frictionalPressure
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{

    return
        Fr*pow(max(alpha - alphaMinFriction, scalar(0)), eta)
       /pow(max(alphaMax - alpha, scalar(5.0e-2)), p);
}


Foam::tmp<Foam::volScalarField> Foam::JohnsonJacksonFrictionalStress::
frictionalPressurePrime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return Fr*
    (
        eta*pow(max(alpha - alphaMinFriction, scalar(0)), eta - 1.0)
       *(alphaMax-alpha) + p*pow(max(alpha - alphaMinFriction, scalar(0)), eta)
    )/pow(max(alphaMax - alpha, scalar(5.0e-2)), p + 1.0);
}


Foam::tmp<Foam::volScalarField> Foam::JohnsonJacksonFrictionalStress::muf
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volTensorField& D,
    const dimensionedScalar& phi
) const
{
    return dimensionedScalar("0.5", dimTime, 0.5)*pf*sin(phi);
}


// ************************************************************************* //
