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

#include "scaleSimilarity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scaleSimilarity, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

scaleSimilarity::scaleSimilarity
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(typeName, U, phi, transport),
    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> scaleSimilarity::k() const
{
    return(0.5*(filter_(magSqr(U())) - magSqr(filter_(U()))));
}


tmp<volScalarField> scaleSimilarity::epsilon() const
{
    volSymmTensorField D = symm(fvc::grad(U()));

    return((filter_(sqr(U())) - sqr(filter_(U()))) && D);
}


tmp<volSymmTensorField> scaleSimilarity::B() const
{
    return(filter_(sqr(U())) - sqr(filter_(U())));
}


tmp<volSymmTensorField> scaleSimilarity::devBeff() const
{
    return dev(B());
}


tmp<fvVectorMatrix> scaleSimilarity::divDevBeff() const
{
    return fvm::Su(fvc::div(devBeff()), U_);
}


void scaleSimilarity::correct(const tmp<volTensorField>&)
{}


bool scaleSimilarity::read()
{
    if (LESModel::read())
    {
        filter_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
