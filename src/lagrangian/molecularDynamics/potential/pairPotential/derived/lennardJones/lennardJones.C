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

#include "lennardJones.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(lennardJones, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    lennardJones,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lennardJones::lennardJones
(
    const word& name,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, pairPotentialProperties),
    lennardJonesCoeffs_(pairPotentialProperties.subDict(typeName + "Coeffs")),
    sigma_(readScalar(lennardJonesCoeffs_.lookup("sigma"))),
    epsilon_(readScalar(lennardJonesCoeffs_.lookup("epsilon")))
{
    setLookupTables();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar lennardJones::unscaledEnergy(const scalar r) const
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;

    return 4.0 * epsilon_*(ir6*(ir6 - 1.0));
}


bool lennardJones::read(const dictionary& pairPotentialProperties)
{
    pairPotential::read(pairPotentialProperties);

    lennardJonesCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    lennardJonesCoeffs_.lookup("sigma") >> sigma_;
    lennardJonesCoeffs_.lookup("epsilon") >> epsilon_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
