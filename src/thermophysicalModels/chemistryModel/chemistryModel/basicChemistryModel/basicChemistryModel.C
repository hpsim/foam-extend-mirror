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

#include "basicChemistryModel.H"
#include "fvMesh.H"
#include "foamTime.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicChemistryModel, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicChemistryModel::basicChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    IOdictionary
    (
        IOobject
        (
            "chemistryProperties",
            mesh.time().constant(),
            obj,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    chemistry_(lookup("chemistry")),
    deltaTChem_
    (
        IOobject
        (
            "deltaTChem",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "deltaTChem0",
            dimTime,
            readScalar(lookup("initialChemicalTimeStep"))
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicChemistryModel::~basicChemistryModel()
{}


// ************************************************************************* //
