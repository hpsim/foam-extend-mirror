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

Description
    Picks up cells with cell centre 'inside' of surface. Requires surface
    to be closed and singly connected.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::noParallel();
    Foam::argList::validArgs.append("surface file");
    Foam::argList::validArgs.append("destination cellSet");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    fileName surfName(args.additionalArgs()[0]);

    fileName setName(args.additionalArgs()[1]);


    // Read surface
    Info<< "Reading surface from " << surfName << endl;
    triSurface surf(surfName);

    // Destination cellSet.
    cellSet insideCells(mesh, setName, IOobject::NO_READ);


    // Construct search engine on surface
    triSurfaceSearch querySurf(surf);

    boolList inside(querySurf.calcInside(mesh.cellCentres()));

    forAll(inside, cellI)
    {
        if (inside[cellI])
        {
            insideCells.insert(cellI);
        }
    }


    Info<< "Selected " << insideCells.size()
        << " cells out of " << mesh.nCells() << endl
        << endl
        << "Writing selected cells to cellSet " << insideCells.name()
        << endl << endl
        << "Use this cellSet e.g. with subsetMesh : " << endl << endl
        << "    subsetMesh <root> <case> " << insideCells.name()
        << endl << endl;

    insideCells.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
