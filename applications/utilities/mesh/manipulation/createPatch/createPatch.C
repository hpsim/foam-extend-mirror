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
    Utility to create patches out of selected boundary faces. Faces come either
    from existing patches or from a faceSet.

    More specifically it:
    - creates new patches (from selected boundary faces). Synchronise faces
      on coupled patches.
    - synchronises points on coupled boundaries
    - remove patches with 0 faces in them

\*---------------------------------------------------------------------------*/

#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "SortableList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"
#include "mapPolyMesh.H"
#include "directTopoChange.H"
#include "polyModifyFace.H"
#include "wordReList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}

// Combine operator to synchronise points. We choose point nearest to origin so
// we can use e.g. great,great,great as null value.
class nearestEqOp
{

public:

    void operator()(vector& x, const vector& y) const
    {
        if (magSqr(y) < magSqr(x))
        {
            x = y;
        }
    }
};


void changePatchID
(
    const polyMesh& mesh,
    const label faceID,
    const label patchID,
    directTopoChange& meshMod
)
{
    const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


// Filter out the empty patches.
void filterPatches(polyMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Patches to keep
    DynamicList<polyPatch*> allPatches(patches.size());

    label nOldPatches = returnReduce(patches.size(), sumOp<label>());

    // Copy old patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Note: reduce possible since non-proc patches guaranteed in same order
        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduce(pp.size(), sumOp<label>()) > 0)
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        allPatches.size(),
                        pp.size(),
                        pp.start()
                    ).ptr()
                );
            }
            else
            {
                Info<< "Removing empty patch " << pp.name() << " at position "
                    << patchI << endl;
            }
        }
    }
    // Copy non-empty processor patches
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            if (pp.size())
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        allPatches.size(),
                        pp.size(),
                        pp.start()
                    ).ptr()
                );
            }
            else
            {
                Info<< "Removing empty processor patch " << pp.name()
                    << " at position " << patchI << endl;
            }
        }
    }

    label nAllPatches = returnReduce(allPatches.size(), sumOp<label>());
    if (nAllPatches != nOldPatches)
    {
        Info<< "Removing zero sizes patches." << endl;
        allPatches.shrink();
        mesh.removeBoundary();
        mesh.addPatches(allPatches);
    }
    else
    {
        Info<< "No patches removed." << endl;
    }
}


// Dump for all patches the current match
void dumpCyclicMatch(const fileName& prefix, const polyMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            label halfSize = cycPatch.size()/2;

            // Dump halves
            {
                OFstream str(prefix+cycPatch.name()+"_half0.obj");
                Pout<< "Dumping " << cycPatch.name()
                    << " half0 faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    static_cast<faceList>
                    (
                        SubList<face>
                        (
                            cycPatch,
                            halfSize
                        )
                    ),
                    cycPatch.points()
                );
            }
            {
                OFstream str(prefix+cycPatch.name()+"_half1.obj");
                Pout<< "Dumping " << cycPatch.name()
                    << " half1 faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    static_cast<faceList>
                    (
                        SubList<face>
                        (
                            cycPatch,
                            halfSize,
                            halfSize
                        )
                    ),
                    cycPatch.points()
                );
            }


            // Lines between corresponding face centres
            OFstream str(prefix+cycPatch.name()+"_match.obj");
            label vertI = 0;

            Pout<< "Dumping cyclic match as lines between face centres to "
                << str.name() << endl;

            for (label faceI = 0; faceI < halfSize; faceI++)
            {
                const point& fc0 = mesh.faceCentres()[cycPatch.start()+faceI];
                meshTools::writeOBJ(str, fc0);
                vertI++;

                label nbrFaceI = halfSize + faceI;
                const point& fc1 =
                    mesh.faceCentres()[cycPatch.start()+nbrFaceI];
                meshTools::writeOBJ(str, fc1);
                vertI++;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
}


void separateList
(
    const vectorField& separation,
    UList<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.

        forAll(field, i)
        {
            field[i] += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] += separation[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "separateList(const vectorField&, UList<vector>&)"
        )   << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


// Synchronise points on both sides of coupled boundaries.
template <class CombineOp>
void syncPoints
(
    const polyMesh& mesh,
    pointField& points,
    const CombineOp& cop,
    const point& nullValue
)
{
    if (points.size() != mesh.nPoints())
    {
        FatalErrorIn
        (
            "syncPoints"
            "(const polyMesh&, pointField&, const CombineOp&, const point&)"
        )   << "Number of values " << points.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if
            (
                isA<processorPolyPatch>(pp)
             && pp.nPoints() > 0
             && refCast<const processorPolyPatch>(pp).master()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Get data per patchPoint in neighbouring point numbers.
                pointField patchInfo(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll(nbrPts, pointI)
                {
                    label nbrPointI = nbrPts[pointI];
                    if (nbrPointI >= 0 && nbrPointI < patchInfo.size())
                    {
                        patchInfo[nbrPointI] = points[meshPts[pointI]];
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
            }
        }


        // Receive and set.

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if
            (
                isA<processorPolyPatch>(pp)
             && pp.nPoints() > 0
             && !refCast<const processorPolyPatch>(pp).master()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                pointField nbrPatchInfo(procPatch.nPoints());
                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }
                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointI)
                {
                    label meshPointI = meshPts[pointI];
                    points[meshPointI] = nbrPatchInfo[pointI];
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();

            pointField half0Values(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point0 = meshPts[e[0]];
                half0Values[i] = points[point0];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
            }
            else if (cycPatch.separated())
            {
                hasTransformation = true;
                const vectorField& v = cycPatch.coupledPolyPatch::separation();
                separateList(v, half0Values);
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point1 = meshPts[e[1]];
                points[point1] = half0Values[i];
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly nessecary.
    //reduce(hasTransformation, orOp<bool>());

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningIn
            (
                "syncPoints"
                "(const polyMesh&, pointField&, const CombineOp&, const point&)"
            )   << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }


        // Values on shared points.
        pointField sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = points[meshPointI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            points[meshPointI] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    argList::validOptions.insert("keepZeroSizedPatches", "");
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();

    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    const bool overwrite = args.optionFound("overwrite");
    const bool keepZeroSizedPatches = args.optionFound("keepZeroSizedPatches");

    Info<< "Reading createPatchDict." << nl << endl;

    IOdictionary dict
    (
        IOobject
        (
            "createPatchDict",
            runTime.system(),
            (
                meshRegionName != polyMesh::defaultRegion
              ? meshRegionName
              : word::null
            ),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    // Whether to synchronise points
    const Switch pointSync(dict.lookup("pointSync"));


    // Change tolerance in controlDict instead.  HJ, 22/Oct/2008

    // Set the matching tolerance so we can read illegal meshes
//     scalar tol = readScalar(dict.lookup("matchTolerance"));
//     Info<< "Using relative tolerance " << tol
//         << " to match up faces and points" << nl << endl;
//      polyPatch::matchTol_ = tol;

#   include "createNamedPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // If running parallel check same patches everywhere
    patches.checkParallelSync(true);


    dumpCyclicMatch("initial_", mesh);

    // Read patch construct info from dictionary
    PtrList<dictionary> patchSources(dict.lookup("patchInfo"));



    // 1. Add all new patches
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (patchSources.size())
    {
        // Old and new patches.
        DynamicList<polyPatch*> allPatches(patches.size()+patchSources.size());

        label startFaceI = mesh.nInternalFaces();

        // Copy old patches.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (!isA<processorPolyPatch>(pp))
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        patchI,
                        pp.size(),
                        startFaceI
                    ).ptr()
                );
                startFaceI += pp.size();
            }
        }

        forAll(patchSources, addedI)
        {
            const dictionary& dict = patchSources[addedI];

            word patchName(dict.lookup("name"));

            label destPatchI = patches.findPatchID(patchName);

            if (destPatchI == -1)
            {
                dictionary patchDict(dict.subDict("dictionary"));

                destPatchI = allPatches.size();

                Info<< "Adding new patch " << patchName
                    << " as patch " << destPatchI
                    << " from " << patchDict << endl;

                patchDict.set("nFaces", 0);
                patchDict.set("startFace", startFaceI);

                // Add an empty patch.
                allPatches.append
                (
                    polyPatch::New
                    (
                        patchName,
                        patchDict,
                        destPatchI,
                        patches
                    ).ptr()
                );
            }
        }

        // Copy old patches.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isA<processorPolyPatch>(pp))
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        patchI,
                        pp.size(),
                        startFaceI
                    ).ptr()
                );
                startFaceI += pp.size();
            }
        }

        allPatches.shrink();
        mesh.removeBoundary();
        mesh.addPatches(allPatches);

        Info<< endl;
    }



    // 2. Repatch faces
    // ~~~~~~~~~~~~~~~~

    directTopoChange meshMod(mesh);


    forAll(patchSources, addedI)
    {
        const dictionary& dict = patchSources[addedI];

        word patchName(dict.lookup("name"));

        label destPatchI = patches.findPatchID(patchName);

        if (destPatchI == -1)
        {
            FatalErrorIn(args.executable()) << "patch " << patchName
                << " not added. Problem." << abort(FatalError);
        }

        word sourceType(dict.lookup("constructFrom"));

        if (sourceType == "patches")
        {
            labelHashSet patchSources
            (
                patches.patchSet(wordReList(dict.lookup("patches")))
            );

            // Repatch faces of the patches.
            forAllConstIter(labelHashSet, patchSources, iter)
            {
                const polyPatch& pp = patches[iter.key()];

                Info<< "Moving faces from patch " << pp.name()
                    << " to patch " << destPatchI << endl;

                forAll(pp, i)
                {
                    changePatchID
                    (
                        mesh,
                        pp.start() + i,
                        destPatchI,
                        meshMod
                    );
                }
            }
        }
        else if (sourceType == "set")
        {
            word setName(dict.lookup("set"));

            faceSet faces(mesh, setName);

            Info<< "Read " << returnReduce(faces.size(), sumOp<label>())
                << " faces from faceSet " << faces.name() << endl;

            // Sort (since faceSet contains faces in arbitrary order)
            labelList faceLabels(faces.toc());

            SortableList<label> patchFaces(faceLabels);

            forAll(patchFaces, i)
            {
                label faceI = patchFaces[i];

                if (mesh.isInternalFace(faceI))
                {
                    FatalErrorIn(args.executable())
                        << "Face " << faceI << " specified in set "
                        << faces.name()
                        << " is not an external face of the mesh." << endl
                        << "This application can only repatch existing boundary"
                        << " faces." << exit(FatalError);
                }

                changePatchID
                (
                    mesh,
                    faceI,
                    destPatchI,
                    meshMod
                );
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Invalid source type " << sourceType << endl
                << "Valid source types are 'patches' 'set'" << exit(FatalError);
        }
    }
    Info<< endl;


    // Change mesh, use inflation to reforce calculation of transformation
    // tensors.
    Info<< "Doing topology modification to order faces." << nl << endl;
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);
    mesh.movePoints(map().preMotionPoints());

    dumpCyclicMatch("coupled_", mesh);

    // Synchronise points.
    if (!pointSync)
    {
        Info<< "Not synchronising points." << nl << endl;
    }
    else
    {
        Info<< "Synchronising points." << nl << endl;

        // This is a bit tricky. Both normal and position might be out and
        // current separation also includes the normal
        // ( separation_ = (nf&(Cr - Cf))*nf ).

        // For processor patches:
        // - disallow multiple separation/transformation. This basically
        //   excludes decomposed cyclics. Use the (probably 0) separation
        //   to align the points.
        // For cyclic patches:
        // - for separated ones use our own recalculated offset vector
        // - for rotational ones use current one.

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            if (pp.size() && isA<coupledPolyPatch>(pp))
            {
                const coupledPolyPatch& cpp =
                    refCast<const coupledPolyPatch>(pp);

                if (cpp.separated())
                {
                    Info<< "On coupled patch " << pp.name()
                        << " separation[0] was "
                        << cpp.separation()[0] << endl;

                    if (isA<cyclicPolyPatch>(pp))
                    {
                        const cyclicPolyPatch& cycpp =
                            refCast<const cyclicPolyPatch>(pp);

                        if (cycpp.transform() == cyclicPolyPatch::TRANSLATIONAL)
                        {
                            Info<< "On cyclic translation patch " << pp.name()
                                << " forcing uniform separation of "
                                << cycpp.separationVector() << endl;
                            const_cast<vectorField&>(cpp.separation()) =
                                pointField(1, cycpp.separationVector());
                        }
                        else
                        {
                            const_cast<vectorField&>(cpp.separation()) =
                                pointField
                                (
                                    1,
                                    pp[pp.size()/2].centre(mesh.points())
                                  - pp[0].centre(mesh.points())
                                );
                        }
                    }
                    else
                    {
                        const_cast<vectorField&>(cpp.separation())
                        .setSize(1);
                    }
                    Info<< "On coupled patch " << pp.name()
                        << " forcing uniform separation of "
                        << cpp.separation() << endl;
                }
                else if (!cpp.parallel())
                {
                    Info<< "On coupled patch " << pp.name()
                        << " forcing uniform rotation of "
                        << cpp.forwardT()[0] << endl;

                    const_cast<tensorField&>
                    (
                        cpp.forwardT()
                    ).setSize(1);
                    const_cast<tensorField&>
                    (
                        cpp.reverseT()
                    ).setSize(1);

                    Info<< "On coupled patch " << pp.name()
                        << " forcing uniform rotation of "
                        << cpp.forwardT() << endl;
                }
            }
        }

        Info<< "Synchronising points." << endl;

        pointField newPoints(mesh.points());

        syncPoints
        (
            mesh,
            newPoints,
            nearestEqOp(),
            point(GREAT, GREAT, GREAT)
        );

        scalarField diff(mag(newPoints-mesh.points()));
        Info<< "Points changed by average:" << gAverage(diff)
            << " max:" << gMax(diff) << nl << endl;

        mesh.movePoints(newPoints);
    }

    // 3. Remove zeros-sized patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!keepZeroSizedPatches)
    {
        Info<< "Removing patches with no faces in them." << nl<< endl;
        filterPatches(mesh);
    }

    dumpCyclicMatch("final_", mesh);


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info<< "Writing repatched mesh to " << runTime.timeName() << nl << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
