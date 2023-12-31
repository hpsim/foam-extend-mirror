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

#include "engineMesh.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::engineMesh, 0);

defineRunTimeSelectionTable(Foam::engineMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::engineMesh::engineMesh(const IOobject& io)
:
    fvMesh(io),
    engineDB_(refCast<const engineTime>(time())),
    pistonIndex_(-1),
    linerIndex_(-1),
    cylinderHeadIndex_(-1),
    deckHeight_("deckHeight", dimLength, GREAT),
    pistonPosition_("deckHeight", dimLength, GREAT)
{
    bool foundPiston = false;
    bool foundLiner = false;
    bool foundCylinderHead = false;

    forAll(boundary(), i)
    {
        if (boundary()[i].name() == "piston")
        {
            pistonIndex_ = i;
            foundPiston = true;
        }
        else if (boundary()[i].name() == "liner")
        {
            linerIndex_ = i;
            foundLiner = true;
        }
        else if (boundary()[i].name() == "cylinderHead")
        {
            cylinderHeadIndex_ = i;
            foundCylinderHead = true;
        }
    }

    reduce(foundPiston, orOp<bool>());
    reduce(foundLiner, orOp<bool>());
    reduce(foundCylinderHead, orOp<bool>());

    if (!foundPiston)
    {
        FatalErrorIn("engineMesh::engineMesh(const IOobject& io)")
            << "cannot find piston patch"
            << exit(FatalError);
    }

    if (!foundLiner)
    {
        FatalErrorIn("engineMesh::engineMesh(const IOobject& io)")
            << "cannot find liner patch"
            << exit(FatalError);
    }

    if (!foundCylinderHead)
    {
        FatalErrorIn("engineMesh::engineMesh(const IOobject& io)")
            << "cannot find cylinderHead patch"
            << exit(FatalError);
    }

    {
        if (pistonIndex_ != -1)
        {
            pistonPosition_.value() =
                max(boundary()[pistonIndex_].patch().localPoints()).z();
        }
        reduce(pistonPosition_.value(), minOp<scalar>());

        if (cylinderHeadIndex_ != -1)
        {
            deckHeight_.value() = min
            (
                boundary()[cylinderHeadIndex_].patch().localPoints()
            ).z();
        }
        reduce(deckHeight_.value(), minOp<scalar>());

        Info<< "deckHeight: " << deckHeight_.value() << nl
            << "piston position: " << pistonPosition_.value() << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::engineMesh::~engineMesh()
{}

// ************************************************************************* //
