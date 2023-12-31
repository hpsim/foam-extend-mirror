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

#include "HrenyaSinclairViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(HrenyaSinclairViscosity, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HrenyaSinclairViscosity,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::HrenyaSinclairViscosity::HrenyaSinclairViscosity
(
    const dictionary& dict
)
:
    viscosityModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    L_(coeffsDict_.lookup("L"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::HrenyaSinclairViscosity::~HrenyaSinclairViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::HrenyaSinclairViscosity::mua
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const volScalarField& g0,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(mathematicalConstant::pi);

    volScalarField lamda =
        scalar(1) + da/(6.0*sqrt(2.0)*(alpha + scalar(1.0e-5)))/L_;

    return rhoa*da*sqrt(Theta)*
    (
        (4.0/5.0)*sqr(alpha)*g0*(1.0 + e)/sqrtPi
      + (1.0/15.0)*sqrtPi*g0*(1.0 + e)*(3.0*e - 1)*sqr(alpha)/(3.0-e)
      + (1.0/6.0)*sqrtPi*alpha*(0.5*lamda + 0.25*(3.0*e - 1.0))
       /(0.5*(3.0 - e)*lamda)
      + (10/96.0)*sqrtPi/((1.0 + e)*0.5*(3.0 - e)*g0*lamda)
    );
}


// ************************************************************************* //
