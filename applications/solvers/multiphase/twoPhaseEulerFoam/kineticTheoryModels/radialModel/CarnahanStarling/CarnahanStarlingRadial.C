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

#include "CarnahanStarlingRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CarnahanStarlingRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        CarnahanStarlingRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CarnahanStarlingRadial::CarnahanStarlingRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CarnahanStarlingRadial::~CarnahanStarlingRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::CarnahanStarlingRadial::g0
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{

    return
        1.0/(1.0 - alpha)
      + 3.0*alpha/(2.0*sqr(1.0 - alpha))
      + sqr(alpha)/(2.0*pow(1.0 - alpha, 3));
}


Foam::tmp<Foam::volScalarField> Foam::CarnahanStarlingRadial::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    return
        - alpha/sqr(1.0 - alpha)
        + (3.0*(1.0 - alpha) + 6.0*sqr(alpha))/(2.0*(1.0 - alpha))
        + (2.0*alpha*(1.0 - alpha) + 3.0*pow(alpha, 3))
         /(2.0*pow(1.0 - alpha, 4));
}


// ************************************************************************* //
