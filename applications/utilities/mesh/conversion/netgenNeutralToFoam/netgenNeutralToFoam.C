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
    Converts neutral file format as written by Netgen v4.4.

    Example:

    9
      1.000000  1.000000  1.000000
      0.000000  1.000000  1.000000
      0.000000  0.000000  1.000000
      1.000000  0.000000  1.000000
      0.000000  1.000000  0.000000
      1.000000  1.000000  0.000000
      1.000000  0.000000  0.000000
      0.000000  0.000000  0.000000
      0.500000  0.500000  0.500000
    12
       1          7        8        9        3
       1          5        9        6        8
       1          5        9        2        1
       1          4        9        7        6
       1          7        8        6        9
       1          4        6        1        9
       1          5        9        8        2
       1          4        1        2        9
       1          1        6        5        9
       1          2        3        4        9
       1          8        9        3        2
       1          4        9        3        7
    12
       1            1        2        4
       1            3        4        2
       2            5        6        8
       2            7        8        6
       3            1        4        6
       3            7        6        4
       5            2        1        5
       5            6        5        1
       5            3        2        8
       5            5        8        2
       6            4        3        7
       6            8        7        3

NOTE:
    - reverse order of boundary faces using geometric test.
      (not very space efficient)
    - order of tet vertices only tested on one file.
    - all patch/cell/vertex numbers offset by one.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "polyPatch.H"
#include "cellModeller.H"
#include "triFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("Neutral file");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName neuFile(args.additionalArgs()[0]);


    IFstream str(neuFile);


    //
    // Read nodes.
    //
    label nNodes(readLabel(str));

    Info<< "nNodes:" << nNodes << endl;


    pointField points(nNodes);

    forAll(points, pointI)
    {
        scalar x,y,z;

        str >> x >> y >> z;

        points[pointI] = point(x, y, z);
    }




    label nTets(readLabel(str));

    Info<< "nTets:" << nTets << endl;

    const cellModel& tet = *(cellModeller::lookup("tet"));

    cellShapeList cells(nTets);

    labelList tetPoints(4);

    forAll(cells, cellI)
    {
        label domain(readLabel(str));

        if (domain != 1)
        {
            WarningIn(args.executable())
                << "Cannot handle multiple domains"
                << nl << "Ignoring domain " << domain << " setting on line "
                << str.lineNumber() << endl;
        }

        tetPoints[1] = readLabel(str) - 1;
        tetPoints[0] = readLabel(str) - 1;
        tetPoints[2] = readLabel(str) - 1;
        tetPoints[3] = readLabel(str) - 1;

        cells[cellI] = cellShape(tet, tetPoints);
    }



    label nFaces(readLabel(str));

    Info<< "nFaces:" << nFaces << endl;

    // Unsorted boundary faces
    faceList boundaryFaces(nFaces);

    // For each boundary faces the Foam patchID
    labelList boundaryPatch(nFaces, -1);

    // Max patch.
    label maxPatch = 0;

    // Boundary faces as three vertices
    HashTable<label, triFace, Hash<triFace> > vertsToBoundary(nFaces);

    forAll(boundaryFaces, faceI)
    {
        label patchI(readLabel(str));

        if (patchI < 0)
        {
            FatalErrorIn(args.executable())
                << "Invalid boundary region number " << patchI
                << " on line " << str.lineNumber()
                << exit(FatalError);
        }


        maxPatch = max(maxPatch, patchI);

        triFace tri(readLabel(str)-1, readLabel(str)-1, readLabel(str)-1);

        // Store boundary face as is for now. Later on reverse it.
        boundaryFaces[faceI].setSize(3);
        boundaryFaces[faceI][0] = tri[0];
        boundaryFaces[faceI][1] = tri[1];
        boundaryFaces[faceI][2] = tri[2];
        boundaryPatch[faceI] = patchI;

        vertsToBoundary.insert(tri, faceI);
    }

    label nPatches = maxPatch + 1;


    // Use hash of points to get owner cell and orient the boundary face.
    // For storage reasons I store the triangles and loop over the cells instead
    // of the other way around (store cells and loop over triangles) though
    // that would be faster.
    forAll(cells, cellI)
    {
        const cellShape& cll = cells[cellI];

        // Get the four (outwards pointing) faces of the cell
        faceList tris(cll.faces());

        forAll(tris, i)
        {
            const face& f = tris[i];

            // Is there any boundary face with same vertices?
            // (uses commutative hash)
            HashTable<label, triFace, Hash<triFace> >::iterator iter =
                vertsToBoundary.find(triFace(f[0], f[1], f[2]));

            if (iter != vertsToBoundary.end())
            {
                label faceI = iter();
                const triFace& tri = iter.key();

                // Determine orientation of tri v.s. cell centre.
                point cc(cll.centre(points));
                point fc(tri.centre(points));
                vector fn(tri.normal(points));

                if (((fc - cc) & fn) < 0)
                {
                    // Boundary face points inwards. Flip.
                    boundaryFaces[faceI] = boundaryFaces[faceI].reverseFace();
                }

                // Done this face so erase from hash
                vertsToBoundary.erase(iter);
            }
        }
    }


    if (vertsToBoundary.size())
    {
        // Didn't find cells connected to boundary faces.
        WarningIn(args.executable())
            << "There are boundary faces without attached cells."
            << "Boundary faces (as triFaces):" << vertsToBoundary.toc()
            << endl;
    }


    // Storage for boundary faces sorted into patches

    faceListList patchFaces(nPatches);

    wordList patchNames(nPatches);

    forAll(patchNames, patchI)
    {
        patchNames[patchI] = word("patch") + name(patchI);
    }

    wordList patchTypes(nPatches, polyPatch::typeName);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = polyPatch::typeName;
    wordList patchPhysicalTypes(nPatches, polyPatch::typeName);

    {
        // Sort boundaryFaces by patch.
        List<DynamicList<face> > allPatchFaces(nPatches);

        forAll(boundaryPatch, faceI)
        {
            label patchI = boundaryPatch[faceI];

            allPatchFaces[patchI].append(boundaryFaces[faceI]);
        }

        Info<< "Patches:" << nl
            << "\tNeutral Boundary\tPatch name\tSize" << nl
            << "\t----------------\t----------\t----" << endl;

        forAll(allPatchFaces, patchI)
        {
            Info<< '\t' << patchI << "\t\t\t"
                << patchNames[patchI] << "\t\t"
                << allPatchFaces[patchI].size() << endl;

            patchFaces[patchI].transfer(allPatchFaces[patchI]);
        }

        Info<< endl;
    }


    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        xferMove(points),
        cells,
        patchFaces,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    Info<< "Writing mesh to " << runTime.constant() << endl << endl;

    mesh.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
