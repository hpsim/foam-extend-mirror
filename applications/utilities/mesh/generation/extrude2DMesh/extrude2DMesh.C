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

#include "extrude2DMesh.H"
#include "polyMesh.H"
#include "directTopoChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(extrude2DMesh, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::extrude2DMesh::extrude2DMesh(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extrude2DMesh::setRefinement
(
    const direction extrudeDir,
    const scalar thickness,
    const label frontPatchI,
    directTopoChange& meshMod
) const
{
    for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        meshMod.addCell
        (
            -1,     //masterPointID,
            -1,     //masterEdgeID,
            -1,     //masterFaceID,
            cellI,  //masterCellID,
            mesh_.cellZones().whichZone(cellI)  //zoneID
        );
    }


    // Generate points
    // ~~~~~~~~~~~~~~~

    forAll(mesh_.points(), pointI)
    {
        meshMod.addPoint
        (
            mesh_.points()[pointI],
            pointI,
            -1,             // zoneID
            true            // inCell
        );
    }

    //Info<< "Adding offsetted points." << nl << endl;
    forAll(mesh_.points(), pointI)
    {
        point newPoint(mesh_.points()[pointI]);
        newPoint[extrudeDir] = thickness;

        meshMod.addPoint
        (
            newPoint,
            pointI,
            -1,             // zoneID
            true            // inCell
        );
    }


    // Generate faces
    // ~~~~~~~~~~~~~~

    const faceList& faces = mesh_.faces();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label zoneID = mesh_.faceZones().whichZone(faceI);
        bool zoneFlip = false;
        if (zoneID != -1)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
        }

        face newFace(4);
        const face& f = faces[faceI];
        newFace[0] = f[0];
        newFace[1] = f[1];
        newFace[2] = f[1]+mesh_.nPoints();
        newFace[3] = f[0]+mesh_.nPoints();

        meshMod.addFace
        (
            newFace,
            mesh_.faceOwner()[faceI],       // own
            mesh_.faceNeighbour()[faceI],   // nei
            -1,                             // masterPointID
            -1,                             // masterEdgeID
            faceI,                          // masterFaceID
            false,                          // flipFaceFlux
            -1,                             // patchID
            zoneID,                         // zoneID
            zoneFlip                        // zoneFlip
        );
    }

    forAll(patches, patchI)
    {
        label startFaceI = patches[patchI].start();
        label endFaceI = startFaceI + patches[patchI].size();

        for (label faceI = startFaceI; faceI < endFaceI; faceI++)
        {
            label zoneID = mesh_.faceZones().whichZone(faceI);
            bool zoneFlip = false;
            if (zoneID != -1)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            face newFace(4);
            const face& f = faces[faceI];
            newFace[0] = f[0];
            newFace[1] = f[1];
            newFace[2] = f[1]+mesh_.nPoints();
            newFace[3] = f[0]+mesh_.nPoints();

            meshMod.addFace
            (
                newFace,
                mesh_.faceOwner()[faceI],       // own
                -1,                             // nei
                -1,                             // masterPointID
                -1,                             // masterEdgeID
                faceI,                          // masterFaceID
                false,                          // flipFaceFlux
                patchI,                         // patchID
                zoneID,                         // zoneID
                zoneFlip                        // zoneFlip
            );
        }
    }


    // Generate front and back faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(mesh_.cells(), cellI)
    {
        const cell& cFaces = mesh_.cells()[cellI];

        // Make a loop out of faces.
        const face& f = faces[cFaces[0]];

        face frontFace(cFaces.size());
        frontFace[0] = f[0];

        label nextPointI = f[1];
        label nextFaceI = cFaces[0];

        for (label i = 1; i < frontFace.size(); i++)
        {
            frontFace[i] = nextPointI;

            // Find face containing pointI
            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];
                if (faceI != nextFaceI)
                {
                    const face& f = faces[faceI];

                    if (f[0] == nextPointI)
                    {
                        nextPointI = f[1];
                        nextFaceI = faceI;
                        break;
                    }
                    else if (f[1] == nextPointI)
                    {
                        nextPointI = f[0];
                        nextFaceI = faceI;
                        break;
                    }
                }
            }
        }


        // Add back face.
        meshMod.addFace
        (
            frontFace.reverseFace(),
            cellI,                          // own
            -1,                             // nei
            -1,                             // masterPointID
            -1,                             // masterEdgeID
            cFaces[0],                      // masterFaceID
            false,                          // flipFaceFlux
            frontPatchI,                    // patchID
            -1,                             // zoneID
            false                           // zoneFlip
        );

        // Offset to create front face.
        forAll(frontFace, fp)
        {
            frontFace[fp] += mesh_.nPoints();
        }
        meshMod.addFace
        (
            frontFace,
            cellI,                          // own
            -1,                             // nei
            -1,                             // masterPointID
            -1,                             // masterEdgeID
            cFaces[0],                      // masterFaceID
            false,                          // flipFaceFlux
            frontPatchI,                    // patchID
            -1,                             // zoneID
            false                           // zoneFlip
        );
    }
}


// ************************************************************************* //
