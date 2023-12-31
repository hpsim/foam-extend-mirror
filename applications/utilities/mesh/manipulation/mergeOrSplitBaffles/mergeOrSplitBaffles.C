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
    mergeOrSplitBaffles

Description
    Detects faces that share points (baffles). Either merge them or
    duplicate the points.

    Notes:
    - can only handle pairwise boundary faces. So three faces using
      the same points is not handled (is illegal mesh anyway)

    - there is no option to only split/merge some baffles.

    - surfaces consisting of duplicate faces can be topologically split
    if the points on the interior of the surface cannot walk to all the
    cells that use them in one go.

    - Parallel operation (where duplicate face is perpendicular to a coupled
    boundary) is supported but not really tested.
    (Note that coupled faces themselves are not seen as duplicate faces)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "syncTools.H"
#include "faceSet.H"
#include "pointSet.H"
#include "meshTools.H"
#include "directTopoChange.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "indirectPrimitivePatch.H"
#include "processorPolyPatch.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void insertDuplicateMerge
(
    const polyMesh& mesh,
    const labelList& duplicates,
    directTopoChange& meshMod
)
{
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const faceZoneMesh& faceZones = mesh.faceZones();

    forAll(duplicates, bFaceI)
    {
        label otherFaceI = duplicates[bFaceI];

        if (otherFaceI != -1 && otherFaceI > bFaceI)
        {
            // Two duplicate faces. Merge.

            label face0 = mesh.nInternalFaces() + bFaceI;
            label face1 = mesh.nInternalFaces() + otherFaceI;

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        own1,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }
    }
}


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces)
{
    // Get all duplicate face labels (in boundaryFaces indices!).
    labelList duplicates = localPointRegion::findDuplicateFaces
    (
        mesh,
        boundaryFaces
    );


    // Check that none are on processor patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(duplicates, bFaceI)
    {
        if (duplicates[bFaceI] != -1)
        {
            label faceI = mesh.nInternalFaces() + bFaceI;
            label patchI = patches.whichPatch(faceI);

            if (isA<processorPolyPatch>(patches[patchI]))
            {
                FatalErrorIn("findBaffles(const polyMesh&, const labelList&)")
                    << "Duplicate face " << faceI
                    << " is on a processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << faceI
                    << " is on patch:" << patches[patchI].name()
                    << abort(FatalError);
            }
        }
    }


    // Write to faceSet for ease of postprocessing.
    {
        faceSet duplicateSet
        (
            mesh,
            "duplicateFaces",
            (mesh.nFaces() - mesh.nInternalFaces())/256
        );

        forAll(duplicates, bFaceI)
        {
            label otherFaceI = duplicates[bFaceI];

            if (otherFaceI != -1 && otherFaceI > bFaceI)
            {
                duplicateSet.insert(mesh.nInternalFaces() + bFaceI);
                duplicateSet.insert(mesh.nInternalFaces() + otherFaceI);
            }
        }

        Pout<< "Writing " << duplicateSet.size()
            << " duplicate faces to faceSet " << duplicateSet.objectPath()
            << nl << endl;
        duplicateSet.write();
    }

    return duplicates;
}




int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    argList::validOptions.insert("split", "");
    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("detectOnly", "");
#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createNamedMesh.H"
    const word oldInstance = mesh.pointsInstance();

    bool split = args.optionFound("split");
    bool overwrite  = args.optionFound("overwrite");
    bool detectOnly = args.optionFound("detectOnly");

    // Collect all boundary faces
    labelList boundaryFaces(mesh.nFaces() - mesh.nInternalFaces());

    forAll (boundaryFaces, i)
    {
        boundaryFaces[i] = i + mesh.nInternalFaces();
    }


    if (detectOnly)
    {
        findBaffles(mesh, boundaryFaces);

        return 0;
    }



    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


    // Mesh change engine
    directTopoChange meshMod(mesh);


    if (split)
    {
        Pout<< "Topologically splitting duplicate surfaces"
            << ", i.e. duplicating points internal to duplicate surfaces."
            << nl << endl;

        // Analyse which points need to be duplicated
        localPointRegion regionSide(mesh);

        // Point duplication engine
        duplicatePoints pointDuplicator(mesh);

        // Insert topo changes
        pointDuplicator.setRefinement(regionSide, meshMod);
    }
    else
    {
        Pout<< "Merging duplicate faces."
            << nl << endl;

        // Get all duplicate face labels (in boundaryFaces indices!).
        labelList duplicates(findBaffles(mesh, boundaryFaces));

        // Merge into internal faces.
        insertDuplicateMerge(mesh, duplicates, meshMod);
    }

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. No inflation.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    Pout<< "Writing mesh to time " << runTime.timeName() << endl;
    mesh.write();

    // Dump duplicated points (if any)
    if (split)
    {
        const labelList& pointMap = map().pointMap();

        labelList nDupPerPoint(map().nOldPoints(), 0);

        pointSet dupPoints(mesh, "duplicatedPoints", 100);

        forAll(pointMap, pointI)
        {
            label oldPointI = pointMap[pointI];

            nDupPerPoint[oldPointI]++;

            if (nDupPerPoint[oldPointI] > 1)
            {
                dupPoints.insert(map().reversePointMap()[oldPointI]);
                dupPoints.insert(pointI);
            }
        }

        Pout<< "Writing " << dupPoints.size()
            << " duplicated points to pointSet "
            << dupPoints.objectPath() << nl << endl;

        dupPoints.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
