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
    Add two surfaces. Does geometric merge on points. Does not check for
    overlapping/intersecting triangles.

    Keeps patches separate by renumbering.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "IFstream.H"
#include "triFace.H"
#include "triFaceList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("Foam surface file");
    argList::validArgs.append("Foam surface file");
    argList::validArgs.append("Foam output file");
    argList::validOptions.insert("points", "pointsFile");
    argList::validOptions.insert("mergeRegions", "");
    argList args(argc, argv);

    fileName inFileName1(args.additionalArgs()[0]);
    fileName inFileName2(args.additionalArgs()[1]);
    fileName outFileName(args.additionalArgs()[2]);

    bool addPoint = args.optionFound("points");
    bool mergeRegions = args.optionFound("mergeRegions");

    if (addPoint)
    {
        Info<< "Reading a surface and adding points from a file"
            << "; merging the points and writing the surface to another file"
            << nl << endl;

        Info<< "Surface  : " << inFileName1<< nl
            << "Points   : " << args.option("points") << nl
            << "Writing  : " << outFileName << nl << endl;
    }
    else
    {
        Info<< "Reading two surfaces"
            << "; merging points and writing the surface to another file"
            << nl << endl;

        if (mergeRegions)
        {
            Info<< "Regions from the two files will get merged" << nl
                << "Do not use this option if you want to keep the regions"
                << " separate" << nl << endl;
        }
        else
        {
            Info<< "Regions from the two files will not get merged" << nl
                << "Regions from " << inFileName2 << " will get offset so"
                << " as not to overlap with the regions in " << inFileName1
                << nl << endl;
        }


        Info<< "Surface1 : " << inFileName1<< nl
            << "Surface2 : " << inFileName2<< nl
            << "Writing  : " << outFileName << nl << endl;
    }

    const triSurface surface1(inFileName1);

    Info<< "Surface1:" << endl;
    surface1.writeStats(Info);
    Info<< endl;

    const pointField& points1 = surface1.points();

    // Final surface
    triSurface combinedSurf;

    if (addPoint)
    {
        IFstream pointsFile(args.option("points"));
        pointField extraPoints(pointsFile);

        Info<< "Additional Points:" << extraPoints.size() << endl;

        vectorField pointsAll(points1);
        label pointI = pointsAll.size();
        pointsAll.setSize(pointsAll.size() + extraPoints.size());

        forAll(extraPoints, i)
        {
            pointsAll[pointI++] = extraPoints[i];
        }

        combinedSurf = triSurface(surface1, surface1.patches(), pointsAll);
    }
    else
    {
        const triSurface surface2(inFileName2);

        Info<< "Surface2:" << endl;
        surface2.writeStats(Info);
        Info<< endl;


        // Make new storage
        List<labelledTri> facesAll(surface1.size() + surface2.size());

        const pointField& points2 = surface2.points();

        vectorField pointsAll(points1.size() + points2.size());


        label pointi = 0;
        // Copy points1 into pointsAll
        forAll(points1, point1i)
        {
            pointsAll[pointi++] = points1[point1i];
        }
        // Add surface2 points
        forAll(points2, point2i)
        {
            pointsAll[pointi++] = points2[point2i];
        }


        label trianglei = 0;

        // Copy triangles1 into trianglesAll
        forAll(surface1, faceI)
        {
            facesAll[trianglei++] = surface1[faceI];
        }
        label nRegions1 = surface1.patches().size();


        if (!mergeRegions)
        {
            Info<< "Surface " << inFileName1 << " has " << nRegions1
                << " regions"
                << nl
                << "All region numbers in " << inFileName2 << " will be offset"
                << " by this amount" << nl << endl;
        }

        // Add (renumbered) surface2 triangles
        forAll(surface2, faceI)
        {
            const labelledTri& tri = surface2[faceI];

            labelledTri& destTri = facesAll[trianglei++];
            destTri[0] = tri[0] + points1.size();
            destTri[1] = tri[1] + points1.size();
            destTri[2] = tri[2] + points1.size();
            if (mergeRegions)
            {
                destTri.region() = tri.region();
            }
            else
            {
                destTri.region() = tri.region() + nRegions1;
            }
        }

        label nRegions2 = surface2.patches().size();

        geometricSurfacePatchList newPatches;

        if (mergeRegions)
        {
            // Overwrite
            newPatches.setSize(max(nRegions1, nRegions2));

            forAll(surface1.patches(), patchI)
            {
                newPatches[patchI] = surface1.patches()[patchI];
            }
            forAll(surface2.patches(), patchI)
            {
                newPatches[patchI] = surface2.patches()[patchI];
            }
        }
        else
        {
            Info<< "Regions from " << inFileName2 << " have been renumbered:"
                << nl
                << "    old\tnew" << nl;

            for (label regionI = 0; regionI < nRegions2; regionI++)
            {
                Info<< "    " << regionI << '\t' << regionI+nRegions1
                    << nl;
            }
            Info<< nl;

            newPatches.setSize(nRegions1 + nRegions2);

            label newPatchI = 0;

            forAll(surface1.patches(), patchI)
            {
                newPatches[newPatchI++] = surface1.patches()[patchI];
            }

            forAll(surface2.patches(), patchI)
            {
                newPatches[newPatchI++] = surface2.patches()[patchI];
            }
        }


        Info<< "New patches:" << nl;
        forAll(newPatches, patchI)
        {
            Info<< "    " << patchI << '\t' << newPatches[patchI].name() << nl;
        }
        Info<< endl;


        // Construct new surface mesh
        combinedSurf = triSurface(facesAll, newPatches, pointsAll);
    }

    // Merge all common points and do some checks
    combinedSurf.cleanup(true);

    Info<< "Merged surface:" << endl;

    combinedSurf.writeStats(Info);

    Info<< endl;

    Info << "Writing : " << outFileName << endl;

    // No need to 'group' while writing since all in correct order anyway.
    combinedSurf.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
