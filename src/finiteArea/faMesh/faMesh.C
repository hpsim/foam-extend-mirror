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

\*---------------------------------------------------------------------------*/

#include "faMesh.H"
#include "faGlobalMeshData.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "IndirectList.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faMeshLduAddressing.H"
#include "wedgeFaPatch.H"
#include "faPatchData.H"
#include "SortableList.H"
#include "controlSwitches.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMesh, 0);
}

Foam::word Foam::faMesh::meshSubDir = "faMesh";

const Foam::debug::optimisationSwitch
Foam::faMesh::quadricsFit_
(
    "quadricsFit",
    0
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::setPrimitiveMeshData()
{
    if (debug)
    {
        Info<< "void faMesh::setPrimitiveMeshData() const : "
            << "Setting primitive data" << endl;
    }

    const indirectPrimitivePatch& bp = patch();

    // Set faMesh edges
    edges_.setSize(bp.nEdges());
    edgeIndex_.setSize(bp.nEdges());

    label edgeI = -1;


    label nIntEdges = bp.nInternalEdges();

    for (label curEdge = 0; curEdge < nIntEdges; curEdge++)
    {
        edges_[++edgeI] = bp.edges()[curEdge];
        edgeIndex_[curEdge] = edgeI;
    }

    forAll (boundary(), patchI)
    {
        const labelList& curP = boundary()[patchI];

        forAll (curP, eI)
        {
            edges_[++edgeI] = bp.edges()[curP[eI]];
            edgeIndex_[curP[eI]] = edgeI;
        }
    }

    nEdges_ = edges_.size();
    nInternalEdges_ = nIntEdges;


    // Set edge owner and neighbour
    edgeOwner_.setSize(nEdges());
    edgeNeighbour_.setSize(nInternalEdges());

    edgeI = -1;

    const labelListList& bpef = bp.edgeFaces();

    for (label curEdge = 0; curEdge < nIntEdges; curEdge++)
    {
        edgeOwner_[++edgeI] = bpef[curEdge][0];

        edgeNeighbour_[edgeI] = bpef[curEdge][1];
    }

    forAll (boundary(), patchI)
    {
        const labelList& curP = boundary()[patchI];

        forAll (curP, eI)
        {
            edgeOwner_[++edgeI] = bpef[curP[eI]][0];
        }
    }

    // Set number of faces
    nFaces_ = bp.size();

    // Set number of points
    nPoints_ = bp.nPoints();
}


void Foam::faMesh::clearGeomNotAreas() const
{
    if (debug)
    {
        Info<< "void faMesh::clearGeomNotAreas() const : "
            << "Clearing geometry" << endl;
    }

    deleteDemandDrivenData(SPtr_);
    deleteDemandDrivenData(patchPtr_);
    deleteDemandDrivenData(patchStartsPtr_);
    deleteDemandDrivenData(LePtr_);
    deleteDemandDrivenData(magLePtr_);
    deleteDemandDrivenData(centresPtr_);
    deleteDemandDrivenData(edgeCentresPtr_);
    deleteDemandDrivenData(faceAreaNormalsPtr_);
    deleteDemandDrivenData(edgeAreaNormalsPtr_);
    deleteDemandDrivenData(pointAreaNormalsPtr_);
    deleteDemandDrivenData(faceCurvaturesPtr_);
    deleteDemandDrivenData(edgeTransformTensorsPtr_);
}


void Foam::faMesh::clearGeom() const
{
    if (debug)
    {
        Info<< "void faMesh::clearGeom() const : "
            << "Clearing geometry" << endl;
    }

    clearGeomNotAreas();
    deleteDemandDrivenData(S0Ptr_);
    deleteDemandDrivenData(S00Ptr_);
    deleteDemandDrivenData(correctPatchPointNormalsPtr_);
}


void Foam::faMesh::clearAddressing() const
{
    if (debug)
    {
        Info<< "void faMesh::clearAddressing() const : "
            << "Clearing addressing" << endl;
    }

    deleteDemandDrivenData(lduPtr_);
}


void Foam::faMesh::clearOut() const
{
    clearGeom();
    clearAddressing();
    deleteDemandDrivenData(globalMeshDataPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh(const polyMesh& pMesh)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, faMesh>(pMesh),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            time().findInstance(meshDir(), "faceLabels"),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            time().findInstance(meshDir(), "faBoundary"),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from IOobject" << endl;
    }

    setPrimitiveMeshData();

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (isFile(pMesh.time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}


// Construct from components without boundary.
Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const labelList& faceLabels
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, faMesh>(pMesh),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        faceLabels
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from components" << endl;
    }
}


// Construct from definition field
Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const fileName& defFile
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, faMesh>(pMesh),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        labelList(0)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from definition file" << endl;
    }

    // Reading faMeshDefinition dictionary
    IOdictionary faMeshDefinition
    (
        IOobject
        (
            defFile,
            mesh().time().constant(),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList polyMeshPatches
    (
        faMeshDefinition.lookup("polyMeshPatches")
    );

    dictionary bndDict = faMeshDefinition.subDict("boundary");

    wordList faPatchNames = bndDict.toc();

    List<faPatchData> faPatches(faPatchNames.size() + 1);

    forAll (faPatchNames, patchI)
    {
        dictionary curPatchDict =
            bndDict.subDict(faPatchNames[patchI]);

        faPatches[patchI].name_ = faPatchNames[patchI];

        faPatches[patchI].type_ =
            word(curPatchDict.lookup("type"));

        faPatches[patchI].ownPolyPatchID_ =
            mesh().boundaryMesh().findPatchID
            (
                word(curPatchDict.lookup("ownerPolyPatch"))
            );

        faPatches[patchI].ngbPolyPatchID_  =
            mesh().boundaryMesh().findPatchID
            (
                word(curPatchDict.lookup("neighbourPolyPatch"))
            );
    }


    // Setting faceLabels list size
    label size = 0;

    labelList patchIDs(polyMeshPatches.size(), -1);

    forAll (polyMeshPatches, patchI)
    {
        patchIDs[patchI] =
            mesh().boundaryMesh().findPatchID(polyMeshPatches[patchI]);

        size += mesh().boundaryMesh()[patchIDs[patchI]].size();
    }

    faceLabels_ = labelList(size, -1);


    // Filling of faceLabels list
    label faceI = -1;

    sort(patchIDs);

    forAll (polyMeshPatches, patchI)
    {
        label start = mesh().boundaryMesh()[patchIDs[patchI]].start();
        label size  = mesh().boundaryMesh()[patchIDs[patchI]].size();

        for (label i = 0; i < size; i++)
        {
            faceLabels_[++faceI] = start + i;
        }
    }


    // Determination of faPatch ID for each boundary edge.
    // Result is in the bndEdgeFaPatchIDs list
    labelList faceCells(faceLabels_.size(), -1);

    forAll (faceCells, faceI)
    {
        label faceID = faceLabels_[faceI];

        faceCells[faceI] = mesh().faceOwner()[faceID];
    }

    labelList meshEdges =
        patch().meshEdges
        (
            mesh().edges(),
            mesh().cellEdges(),
            faceCells
        );

    const labelListList& edgeFaces = mesh().edgeFaces();

    const label nTotalEdges = patch().nEdges();
    const label nInternalEdges = patch().nInternalEdges();

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
    {
        label curMeshEdge = meshEdges[edgeI];

        labelList curEdgePatchIDs(2, label(-1));

        label patchI = -1;

        forAll (edgeFaces[curMeshEdge], faceI)
        {
            label curFace = edgeFaces[curMeshEdge][faceI];

            label curPatchID = mesh().boundaryMesh().whichPatch(curFace);

            if (curPatchID != -1)
            {
                curEdgePatchIDs[++patchI] = curPatchID;
            }
        }

        for (label pI = 0; pI < faPatches.size() - 1; pI++)
        {
            if
            (
                (
                    curEdgePatchIDs[0] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[1] == faPatches[pI].ngbPolyPatchID_
                )
             ||
                (
                    curEdgePatchIDs[1] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[0] == faPatches[pI].ngbPolyPatchID_
                )
            )
            {
                bndEdgeFaPatchIDs[edgeI - nInternalEdges] = pI;
                break;
            }
        }
    }


    // Set edgeLabels for each faPatch
    for (label pI = 0; pI < (faPatches.size() - 1); pI++)
    {
        SLList<label> tmpList;

        forAll (bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI] == pI)
            {
                tmpList.append(nInternalEdges + eI);
            }
        }

        faPatches[pI].edgeLabels_ = tmpList;
    }

    // Check for undefined edges
    SLList<label> tmpList;

    forAll (bndEdgeFaPatchIDs, eI)
    {
        if (bndEdgeFaPatchIDs[eI] == -1)
        {
            tmpList.append(nInternalEdges + eI);
        }
    }

    if (tmpList.size() > 0)
    {
        // Check for processor edges
        labelList allUndefEdges = tmpList;
        labelList ngbPolyPatch(allUndefEdges.size(), -1);
        forAll (ngbPolyPatch, edgeI)
        {
            label curEdge = allUndefEdges[edgeI];

            label curPMeshEdge = meshEdges[curEdge];

            forAll (edgeFaces[curPMeshEdge], faceI)
            {
                label curFace = edgeFaces[curPMeshEdge][faceI];

                if (findIndex(faceLabels_, curFace) == -1)
                {
                    label polyPatchID =
                        pMesh.boundaryMesh().whichPatch(curFace);

                    if (polyPatchID != -1)
                    {
                        ngbPolyPatch[edgeI] = polyPatchID;
                    }
                }
            }
        }

        // Count ngb processorPolyPatch-es
        labelHashSet processorPatchSet;
        forAll (ngbPolyPatch, edgeI)
        {
            if (ngbPolyPatch[edgeI] != -1)
            {
                if
                (
                    isA<processorPolyPatch>
                    (
                        pMesh.boundaryMesh()[ngbPolyPatch[edgeI]]
                    )
                )
                {
                    if (!processorPatchSet.found(ngbPolyPatch[edgeI]))
                    {
                        processorPatchSet.insert(ngbPolyPatch[edgeI]);
                    }
                }
            }
        }
        labelList processorPatches(processorPatchSet.toc());
        faPatches.setSize(faPatches.size() + processorPatches.size());

        for (label i = 0; i < processorPatches.size(); i++)
        {
            SLList<label> tmpLst;

            forAll (ngbPolyPatch, eI)
            {
                if (ngbPolyPatch[eI] == processorPatches[i])
                {
                    tmpLst.append(allUndefEdges[eI]);
                }
            }

            faPatches[faPatchNames.size() + i].edgeLabels_ = tmpLst;

            faPatches[faPatchNames.size() + i].name_ =
                pMesh.boundaryMesh()[processorPatches[i]].name();

            faPatches[faPatchNames.size() + i].type_ =
                processorFaPatch::typeName;

            faPatches[faPatchNames.size() + i].ngbPolyPatchID_ =
                processorPatches[i];
        }

        // Remaining undefined edges
        SLList<label> undefEdges;
        forAll (ngbPolyPatch, eI)
        {
            if (ngbPolyPatch[eI] == -1)
            {
                undefEdges.append(allUndefEdges[eI]);
            }
            else if
            (
               !isA<processorPolyPatch>
                (
                    pMesh.boundaryMesh()[ngbPolyPatch[eI]]
                )
            )
            {
                undefEdges.append(allUndefEdges[eI]);
            }
        }

        if (undefEdges.size() > 0)
        {
            label pI = faPatches.size()-1;

            faPatches[pI].name_ = "undefined";
            faPatches[pI].type_ = "patch";
            faPatches[pI].edgeLabels_ = undefEdges;
        }
        else
        {
            faPatches.setSize(faPatches.size() - 1);
        }
    }
    else
    {
        faPatches.setSize(faPatches.size() - 1);
    }


    // Reorder processorFaPatch using
    // ordering of ngb processorPolyPatch
    forAll (faPatches, patchI)
    {
        if (faPatches[patchI].type_ == processorFaPatch::typeName)
        {
            labelList ngbFaces(faPatches[patchI].edgeLabels_.size(), -1);

            forAll (ngbFaces, edgeI)
            {
                label curEdge = faPatches[patchI].edgeLabels_[edgeI];

                label curPMeshEdge = meshEdges[curEdge];

                forAll (edgeFaces[curPMeshEdge], faceI)
                {
                    label curFace = edgeFaces[curPMeshEdge][faceI];

                    label curPatchID =
                        pMesh.boundaryMesh().whichPatch(curFace);

                    if (curPatchID == faPatches[patchI].ngbPolyPatchID_)
                    {
                        ngbFaces[edgeI] = curFace;
                    }
                }
            }

            SortableList<label> sortedNgbFaces(ngbFaces);
            labelList reorderedEdgeLabels(ngbFaces.size(), -1);
            for (label i = 0; i < reorderedEdgeLabels.size(); i++)
            {
                reorderedEdgeLabels[i] =
                    faPatches[patchI].edgeLabels_
                    [
                        sortedNgbFaces.indices()[i]
                    ];
            }

            faPatches[patchI].edgeLabels_ = reorderedEdgeLabels;
        }
    }


    // Add good patches to faMesh
    SLList<faPatch*> faPatchLst;

    for (label pI = 0; pI < faPatches.size(); pI++)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        faPatches[pI].dict_.add
        (
            "ngbPolyPatchIndex",
            faPatches[pI].ngbPolyPatchID_
        );

        if (faPatches[pI].type_ == processorFaPatch::typeName)
        {
            if (faPatches[pI].ngbPolyPatchID_ == -1)
            {
                FatalErrorIn
                (
                    "void faMesh::faMesh(const polyMesh&, const fileName&)"
                )
                    << "ngbPolyPatch is not defined for processorFaPatch: "
                        << faPatches[pI].name_
                        << abort(FatalError);
            }

            const processorPolyPatch& procPolyPatch =
                refCast<const processorPolyPatch>
                (
                    pMesh.boundaryMesh()[faPatches[pI].ngbPolyPatchID_]
                );

            faPatches[pI].dict_.add("myProcNo", procPolyPatch.myProcNo());
            faPatches[pI].dict_.add
            (
                "neighbProcNo",
                procPolyPatch.neighbProcNo()
            );
        }

        faPatchLst.append
        (
            faPatch::New
            (
                faPatches[pI].name_,
                faPatches[pI].dict_,
                pI,
                boundary()
            ).ptr()
        );
    }

    addFaPatches(List<faPatch*>(faPatchLst));

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (isFile(mesh().time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}


// Construct from polyPatch
Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const label polyPatchID
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, faMesh>(pMesh),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        labelList(pMesh.boundaryMesh()[polyPatchID].size(), -1)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from polyPatch" << endl;
    }

    // Set faceLabels
    forAll (faceLabels_, faceI)
    {
        faceLabels_[faceI] =
            mesh().boundaryMesh()[polyPatchID].start() + faceI;
    }

    // Add one faPatch
    const indirectPrimitivePatch& bp = patch();

    const label nTotalEdges = bp.nEdges();

    const label nInternalEdges = bp.nInternalEdges();

    labelList edgeLabels(nTotalEdges - nInternalEdges, -1);

    forAll (edgeLabels, edgeI)
    {
        edgeLabels[edgeI] = nInternalEdges + edgeI;
    }

    dictionary patchDict;

    patchDict.add("type", "patch");
    patchDict.add("edgeLabels", edgeLabels);
    patchDict.add("ngbPolyPatchIndex", -1);

    List<faPatch*> faPatchLst(1);

    faPatchLst[0] =
        faPatch::New("default", patchDict, 0, boundary()).ptr();

    addFaPatches(faPatchLst);

    setPrimitiveMeshData();

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faMesh::~faMesh()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::faMesh::meshDir() const
{
    return mesh().dbDir()/meshSubDir;
}


const Foam::Time& Foam::faMesh::time() const
{
    return mesh().time();
}


const Foam::fileName& Foam::faMesh::pointsInstance() const
{
    return mesh().pointsInstance();
}


const Foam::fileName& Foam::faMesh::facesInstance() const
{
    return mesh().facesInstance();
}


const Foam::indirectPrimitivePatch& Foam::faMesh::patch() const
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().allFaces(),
                faceLabels_
            ),
            mesh().allPoints()
        );
    }

    return *patchPtr_;
}


Foam::indirectPrimitivePatch& Foam::faMesh::patch()
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().allFaces(),
                faceLabels_
            ),
            mesh().allPoints()
        );
    }

    return *patchPtr_;
}


const Foam::pointField& Foam::faMesh::points() const
{
    return patch().localPoints();
}


const Foam::edgeList& Foam::faMesh::edges() const
{
    return edges_;
}


const Foam::faceList& Foam::faMesh::faces() const
{
    return patch().localFaces();
}


void Foam::faMesh::addFaPatches(const List<faPatch*>& p)
{
    if (debug)
    {
        Info<< "void faMesh::addFaPatches(const List<faPatch*>& p) : "
            << "Adding patches to faMesh" << endl;
    }

    if (boundary().size() > 0)
    {
        FatalErrorIn("void faMesh::addPatches(const List<faPatch*>& p)")
            << "boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    forAll (p, patchI)
    {
        boundary_.set(patchI, p[patchI]);
    }

    setPrimitiveMeshData();

    boundary_.checkDefinition();
}


int Foam::faMesh::comm() const
{
    return comm_;
}


int& Foam::faMesh::comm()
{
    return comm_;
}


const Foam::objectRegistry& Foam::faMesh::thisDb() const
{
    return mesh().thisDb();
}


const Foam::faBoundaryMesh& Foam::faMesh::boundary() const
{
    return boundary_;
}


const Foam::labelList& Foam::faMesh::patchStarts() const
{
    if (!patchStartsPtr_)
    {
        calcPatchStarts();
    }

    return *patchStartsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::Le() const
{
    if (!LePtr_)
    {
        calcLe();
    }

    return *LePtr_;
}


const Foam::edgeScalarField& Foam::faMesh::magLe() const
{
    if (!magLePtr_)
    {
        calcMagLe();
    }

    return *magLePtr_;
}


const Foam::areaVectorField& Foam::faMesh::areaCentres() const
{
    if (!centresPtr_)
    {
        calcAreaCentres();
    }

    return *centresPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeCentres() const
{
    if (!edgeCentresPtr_)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S() const
{
    if (!SPtr_)
    {
        calcS();
    }

    return *SPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S0() const
{
    if (!S0Ptr_)
    {
        FatalErrorIn("faMesh::S0() const")
            << "S0 is not available"
            << abort(FatalError);
    }

    return *S0Ptr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S00() const
{
    if (!S00Ptr_)
    {
        S00Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S00",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            S0()
        );

        S0Ptr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    return *S00Ptr_;
}


const Foam::areaVectorField& Foam::faMesh::faceAreaNormals() const
{
    if (!faceAreaNormalsPtr_)
    {
        calcFaceAreaNormals();
    }

    return *faceAreaNormalsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeAreaNormals() const
{
    if (!edgeAreaNormalsPtr_)
    {
        calcEdgeAreaNormals();
    }

    return *edgeAreaNormalsPtr_;
}


const Foam::vectorField& Foam::faMesh::pointAreaNormals() const
{
    if (!pointAreaNormalsPtr_)
    {
        calcPointAreaNormals();

        if (quadricsFit_() > 0)
        {
            calcPointAreaNormalsByQuadricsFit();
        }
    }

    return *pointAreaNormalsPtr_;
}


const Foam::areaScalarField& Foam::faMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        calcFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


const Foam::FieldField<Foam::Field, Foam::tensor>&
Foam::faMesh::edgeTransformTensors() const
{
    if (!edgeTransformTensorsPtr_)
    {
        calcEdgeTransformTensors();
    }

    return *edgeTransformTensorsPtr_;
}


// Return global mesh data
const Foam::faGlobalMeshData& Foam::faMesh::globalData() const
{
    if (!globalMeshDataPtr_)
    {
        globalMeshDataPtr_ = new faGlobalMeshData(*this);
    }

    return *globalMeshDataPtr_;
}


const Foam::lduAddressing& Foam::faMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        calcLduAddressing();
    }

    return *lduPtr_;
}


bool Foam::faMesh::movePoints() const
{
    if (debug)
    {
        InfoIn("bool faMesh::movePoints() const")
            << "Moving points" << endl;
    }

    // Grab point motion from polyMesh
    const vectorField& newPoints = mesh().allPoints();

    // Grab old time areas if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (S00Ptr_ && S0Ptr_)
        {
            *S00Ptr_ = *S0Ptr_;
        }

        if (S0Ptr_)
        {
            *S0Ptr_ = S();
        }
        else
        {
            if (debug)
            {
                InfoIn("bool faMesh::movePoints() const")
                    << "Creating old face areas." << endl;
            }

            S0Ptr_ = new DimensionedField<scalar, areaMesh>
            (
                IOobject
                (
                    "S0",
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                S()
            );
        }

        curTimeIndex_ = time().timeIndex();
    }

    clearGeomNotAreas();

    // To satisfy the motion interface for MeshObject, const cast is needed
    // HJ, 5/Aug/2011
    if (patchPtr_)
    {
        patchPtr_->movePoints(newPoints);
    }

    // Move boundary points
    const_cast<faBoundaryMesh&>(boundary_).movePoints(newPoints);

    // Move interpolation
    const edgeInterpolation& cei = *this;
    const_cast<edgeInterpolation&>(cei).edgeInterpolation::movePoints();

    // Fluxes were dummy?  HJ, 28/Jul/2011

    return true;
}


bool Foam::faMesh::correctPatchPointNormals(const label patchID) const
{
    if ((patchID < 0) || (patchID >= boundary().size()))
    {
        FatalErrorIn
        (
            "bool correctPatchPointNormals(const label patchID) const"
        )   << "patchID is not in the valid range"
            << abort(FatalError);
    }

    if (correctPatchPointNormalsPtr_)
    {
        return (*correctPatchPointNormalsPtr_)[patchID];
    }

    return false;
}


//- Set patch point normals corrections
Foam::boolList& Foam::faMesh::correctPatchPointNormals() const
{
    if (!correctPatchPointNormalsPtr_)
    {
        correctPatchPointNormalsPtr_ =
            new boolList(boundary().size(), false);
    }

    return *correctPatchPointNormalsPtr_;
}


bool Foam::faMesh::write() const
{
    faceLabels_.write();
    boundary_.write();

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::faMesh::operator!=(const faMesh& m) const
{
    return &m != this;
}


bool Foam::faMesh::operator==(const faMesh& m) const
{
    return &m == this;
}


// ************************************************************************* //
