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

#include "mixedSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mixedSmagorinsky, 0);
addToRunTimeSelectionTable(LESModel, mixedSmagorinsky, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixedSmagorinsky::mixedSmagorinsky
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    scaleSimilarity(U, phi, transport),
    Smagorinsky(U, phi, transport)
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> mixedSmagorinsky::k() const
{
    return
    (
        scaleSimilarity::k()
      + Smagorinsky::k()
    );
}


tmp<volScalarField> mixedSmagorinsky::epsilon() const
{
    return
    (
        scaleSimilarity::epsilon()
      + Smagorinsky::epsilon()
    );
}


tmp<volSymmTensorField> mixedSmagorinsky::B() const
{
    return
    (
        scaleSimilarity::B()
      + Smagorinsky::B()
    );
}


tmp<volSymmTensorField> mixedSmagorinsky::devBeff() const
{
    return
    (
        scaleSimilarity::devBeff()
      + Smagorinsky::devBeff()
    );
}


tmp<fvVectorMatrix> mixedSmagorinsky::divDevBeff() const
{
    return
    (
        scaleSimilarity::divDevBeff()
      + Smagorinsky::divDevBeff()
    );
}


void mixedSmagorinsky::correct(const tmp<volTensorField>& gradU)
{
    scaleSimilarity::correct(gradU);
    Smagorinsky::correct(gradU);
}


bool mixedSmagorinsky::read()
{
    if (LESModel::read())
    {
        scaleSimilarity::read();
        Smagorinsky::read();

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
