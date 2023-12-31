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

Application
    surfaceAutoPatch

Description
    Patches surface according to feature angle. Like autoPatch.

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "argList.H"
#include "surfaceFeatures.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("included angle [0..180]");
    argList args(argc, argv);

    fileName inFileName(args.additionalArgs()[0]);
    fileName outFileName(args.additionalArgs()[1]);
    scalar includedAngle(readScalar(IStringStream(args.additionalArgs()[2])()));

    Pout<< "Surface        : " << inFileName << nl
        << endl;


    // Read
    // ~~~~

    Info << "Reading : " << inFileName << endl;
    triSurface surf(inFileName);

    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;



    // Construct features from surface&featureangle
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Constructing feature set from included angle " << includedAngle
        << endl;

    surfaceFeatures set(surf, includedAngle);

    Pout<< nl
        << "Feature set:" << nl
        << "    feature points : " << set.featurePoints().size() << nl
        << "    feature edges  : " << set.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << set.nRegionEdges() << nl
        << "        external edges : " << set.nExternalEdges() << nl
        << "        internal edges : " << set.nInternalEdges() << nl
        << endl;

    // Get per-edge status.
    boolList borderEdge(surf.nEdges(), false);
    forAll(set.featureEdges(), i)
    {
        borderEdge[set.featureEdges()[i]] = true;
    }

    labelList faceRegion(surf.size());
    label nRegions = surf.markZones(borderEdge, faceRegion);

    // Reregion triangles.
    forAll(surf, i)
    {
        surf[i].region() = faceRegion[i];
    }

    // Create some patches
    surf.patches().setSize(nRegions);

    forAll(surf.patches(), patchI)
    {
        surf.patches()[patchI].name() = "patch" + Foam::name(patchI);
        surf.patches()[patchI].geometricType() = "empty";
    }


    Info << "Writing : " << outFileName << endl;
    surf.write(outFileName, true);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
