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

#include "pairPatchAgglomeration.H"
#include "meshTools.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairPatchAgglomeration::compactLevels(const label nCreatedLevels)
{
    nFaces_.setSize(nCreatedLevels);
    restrictAddressing_.setSize(nCreatedLevels);
    patchLevels_.setSize(nCreatedLevels);
}


bool Foam::pairPatchAgglomeration::continueAgglomerating
(
    const label nCoarseFaces
)
{
    // Check the need for further agglomeration on all processors
    bool contAgg = nCoarseFaces >= nFacesInCoarsestLevel_;
    reduce(contAgg, andOp<bool>());
    return contAgg;
}


void Foam::pairPatchAgglomeration::setBasedEdgeWeights()
{
    const standAlonePatch& coarsePatch = patchLevels_[0];
    forAll(coarsePatch.edges(), i)
    {
        if (coarsePatch.isInternalEdge(i))
        {
            scalar edgeLength =
                coarsePatch.edges()[i].mag(coarsePatch.localPoints());

            const labelList& eFaces = coarsePatch.edgeFaces()[i];

            if (eFaces.size() == 2)
            {
                // Calculate dot product of face normals
                scalar cosI =
                (
                    coarsePatch.faceNormals()[eFaces[0]]
                  & coarsePatch.faceNormals()[eFaces[1]]
                );

                const edge edgeCommon = edge(eFaces[0], eFaces[1]);

                if (facePairWeight_.found(edgeCommon))
                {
                    facePairWeight_[edgeCommon] += edgeLength;
                }
                else
                {
                    facePairWeight_.insert(edgeCommon, edgeLength);
                }

                if (cosI < Foam::cos(degToRad(featureAngle_)))
                {
                    facePairWeight_[edgeCommon] = -1.0;
                }
            }
            else if (eFaces.size() == 3)
            {
                facePairWeight_.insert(edge(eFaces[1], eFaces[0]), -1.0);
                facePairWeight_.insert(edge(eFaces[2], eFaces[0]), -1.0);
                facePairWeight_.insert(edge(eFaces[1], eFaces[2]), -1.0);
            }
        }
    }
}


void Foam::pairPatchAgglomeration::setEdgeWeights
(
    const label fineLevelIndex
)
{

    const standAlonePatch& coarsePatch = patchLevels_[fineLevelIndex];

    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];
    const label nCoarseI =  max(fineToCoarse) + 1;
    labelListList coarseToFine(invertOneToMany(nCoarseI, fineToCoarse));

    HashSet<edge, Hash<edge> > fineFeaturedFaces(coarsePatch.nEdges()/10);

    // Map fine faces with featured edge into coarse faces
    forAllConstIter(EdgeMap<scalar>, facePairWeight_, iter)
    {
        if (iter() == -1.0)
        {
            const edge e = iter.key();
            const edge edgeFeatured
            (
                fineToCoarse[e[0]],
                fineToCoarse[e[1]]
            );
            fineFeaturedFaces.insert(edgeFeatured);
        }
    }

    // Clean old weitghs
    facePairWeight_.clear();
    facePairWeight_.resize(coarsePatch.nEdges());

    forAll(coarsePatch.edges(), i)
    {
        if (coarsePatch.isInternalEdge(i))
        {
            scalar edgeLength =
                coarsePatch.edges()[i].mag(coarsePatch.localPoints());

            const labelList& eFaces = coarsePatch.edgeFaces()[i];

            if (eFaces.size() == 2)
            {
                const edge edgeCommon = edge(eFaces[0], eFaces[1]);
                if (facePairWeight_.found(edgeCommon))
                {
                    facePairWeight_[edgeCommon] += edgeLength;
                }
                else
                {
                    facePairWeight_.insert(edgeCommon, edgeLength);
                }
                // If the fine 'pair' faces was featured edge so it is
                // the coarse 'pair'
                if (fineFeaturedFaces.found(edgeCommon))
                {
                    facePairWeight_[edgeCommon] = -1.0;
                }
            }
            else if (eFaces.size() == 3)
            {
                facePairWeight_.insert(edge(eFaces[1], eFaces[0]), -1.0);
                facePairWeight_.insert(edge(eFaces[2], eFaces[0]), -1.0);
                facePairWeight_.insert(edge(eFaces[1], eFaces[2]), -1.0);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPatchAgglomeration::pairPatchAgglomeration
(
    const polyPatch& patch,
    const dictionary& controlDict,
    const bool additionalWeights
)
:
    mergeLevels_
    (
        controlDict.lookupOrDefault<label>("mergeLevels", 2)
    ),
    maxLevels_(50),
    nFacesInCoarsestLevel_
    (
        readLabel(controlDict.lookup("nFacesInCoarsestLevel"))
    ),
    featureAngle_
    (
        controlDict.lookupOrDefault<scalar>("featureAngle", 0)
    ),
    nFaces_(maxLevels_),
    restrictAddressing_(maxLevels_),
    restrictTopBottomAddressing_(identity(patch.size())),
    patchLevels_(maxLevels_),
    facePairWeight_(patch.size())
{
    // Set base fine patch
    patchLevels_.set
    (
        0,
        new standAlonePatch
        (
            patch.localFaces(),
            patch.localPoints()
        )
    );

    // Set number of faces for the base patch
    nFaces_[0] = patch.size();

    // Set edge weights for level 0
    setBasedEdgeWeights();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairPatchAgglomeration::~pairPatchAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::standAlonePatch& Foam::pairPatchAgglomeration::patchLevel
(
    const label i
) const
{
    return patchLevels_[i];
}


void Foam::pairPatchAgglomeration::mapBaseToTopAgglom
(
    const label fineLevelIndex
)
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];
    forAll (restrictTopBottomAddressing_, i)
    {
        restrictTopBottomAddressing_[i] =
            fineToCoarse[restrictTopBottomAddressing_[i]];
    }
}


bool Foam::pairPatchAgglomeration::agglomeratePatch
(
    const standAlonePatch& patch,
    const labelList& fineToCoarse,
    const label fineLevelIndex
)
{
    if (min(fineToCoarse) == -1)
    {
        FatalErrorIn
        (
            "pairPatchAgglomeration::agglomeratePatch"
            "("
                "const standAlonePatch&, "
                "const labelList&, "
                "const label"
            ")"
        )   << "min(fineToCoarse) == -1" << exit(FatalError);
    }

    if (fineToCoarse.size() != patch.size())
    {
        FatalErrorIn
        (
            "pairPatchAgglomeration::agglomeratePatch"
            "("
                "const standAlonePatch&, "
                "const labelList&, "
                "const label"
            ")"
        )   << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << fineToCoarse.size()
            << " nEqns: " << patch.size()
            << abort(FatalError);
    }

    const label nCoarseI =  max(fineToCoarse)+1;
    List<face> patchFaces(nCoarseI);

    // Patch faces per agglomeration
    labelListList coarseToFine(invertOneToMany(nCoarseI, fineToCoarse));

    for (label coarseI = 0; coarseI < nCoarseI; coarseI++)
    {
        const labelList& fineFaces = coarseToFine[coarseI];

        // Construct single face
        indirectPrimitivePatch upp
        (
            IndirectList<face>(patch, fineFaces),
            patch.points()
        );

        if (upp.edgeLoops().size() != 1)
        {
            if (fineFaces.size() == 2)
            {
                const edge e(fineFaces[0], fineFaces[1]);
                facePairWeight_[e] = -1.0;
            }
            else if (fineFaces.size() == 3)
            {
                const edge e(fineFaces[0], fineFaces[1]);
                const edge e1(fineFaces[0], fineFaces[2]);
                const edge e2(fineFaces[2], fineFaces[1]);
                facePairWeight_[e] = -1.0;
                facePairWeight_[e1] = -1.0;
                facePairWeight_[e2] = -1.0;
            }

            return false;
        }

        patchFaces[coarseI] = face
        (
            renumber
            (
                upp.meshPoints(),
                upp.edgeLoops()[0]
            )
        );
    }

    patchLevels_.set
    (
        fineLevelIndex,
        new standAlonePatch
        (
            SubList<face>(patchFaces, nCoarseI, 0),
            patch.points()
        )
    );
    return true;
}


void Foam::pairPatchAgglomeration:: agglomerate()
{
    label nPairLevels = 0;
    label nCreatedLevels = 1; //0 level is the base patch

    while (nCreatedLevels < maxLevels_)
    {
        label nCoarseCells = -1;

        const standAlonePatch& patch = patchLevels_[nCreatedLevels - 1];
        tmp<labelField> finalAgglomPtr(new labelField(patch.size()));

        bool agglomOK = false;
        while (!agglomOK && patch.size())
        {
            finalAgglomPtr = agglomerateOneLevel
            (
                nCoarseCells,
                patch
            );

            agglomOK = agglomeratePatch
            (
                patch,
                finalAgglomPtr,
                nCreatedLevels
            );
        }

        nFaces_[nCreatedLevels] = nCoarseCells;
        restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        mapBaseToTopAgglom(nCreatedLevels);

        if (!continueAgglomerating(nCoarseCells))
        {
            break;
        }

        setEdgeWeights(nCreatedLevels);

        if (nPairLevels % mergeLevels_)
        {
            combineLevels(nCreatedLevels);
        }
        else
        {
            nCreatedLevels++;
        }

        nPairLevels++;
    }
}


Foam::tmp<Foam::labelField> Foam::pairPatchAgglomeration::agglomerateOneLevel
(
    label& nCoarseCells,
    const standAlonePatch& patch
)
{
    const label nFineFaces = patch.size();

    tmp<labelField> tcoarseCellMap(new labelField(nFineFaces, -1));
    labelField& coarseCellMap = tcoarseCellMap();

    const labelListList& faceFaces = patch.faceFaces();

    nCoarseCells = 0;

    forAll(faceFaces, facei)
    {
        const labelList& fFaces = faceFaces[facei];

        if (coarseCellMap[facei] < 0)
        {
            label matchFaceNo = -1;
            label matchFaceNeibNo = -1;
            scalar maxFaceWeight = -GREAT;

            // check faces to find ungrouped neighbour with largest face weight
            forAll(fFaces, i)
            {
                label faceNeig = fFaces[i];
                const edge edgeCommon = edge(facei, faceNeig);
                if
                (
                    facePairWeight_[edgeCommon] > maxFaceWeight
                 && coarseCellMap[faceNeig] < 0
                 && facePairWeight_[edgeCommon] != -1.0
                )
                {
                    // Match found. Pick up all the necessary data
                    matchFaceNo = facei;
                    matchFaceNeibNo = faceNeig;
                    maxFaceWeight = facePairWeight_[edgeCommon];
                }
            }

            if (matchFaceNo >= 0)
            {
                // Make a new group
                coarseCellMap[matchFaceNo] = nCoarseCells;
                coarseCellMap[matchFaceNeibNo] = nCoarseCells;
                nCoarseCells++;
            }
            else
            {
                // No match. Find the best neighbouring cluster and
                // put the cell there
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = -GREAT;

                forAll(fFaces, i)
                {
                    label faceNeig = fFaces[i];
                    const edge edgeCommon = edge(facei, faceNeig);
                    if
                    (
                        facePairWeight_[edgeCommon] > clusterMaxFaceCoeff
                        && facePairWeight_[edgeCommon] != -1.0
                        && coarseCellMap[faceNeig] > 0
                    )
                    {
                        clusterMatchFaceNo = faceNeig;
                        clusterMaxFaceCoeff = facePairWeight_[edgeCommon];
                    }
                }

                if (clusterMatchFaceNo >= 0)
                {
                    // Add the cell to the best cluster
                    coarseCellMap[facei] = coarseCellMap[clusterMatchFaceNo];
                }
                else
                {
                    // if not create single-cell "clusters" for each
                    coarseCellMap[facei] = nCoarseCells;
                    nCoarseCells ++;
                }
            }
        }
    }

    // Check that all faces are part of clusters,

    for (label facei = 0; facei < nFineFaces; facei++)
    {
        if (coarseCellMap[facei] < 0)
        {
            FatalErrorIn
            (
                "pairPatchAgglomeration::agglomerateOneLevel"
                "(label&, const standAlonePatch&) "
            ) << " face " << facei
            << " is not part of a cluster"
            << exit(FatalError);
        }
    }

    return tcoarseCellMap;
}


void Foam::pairPatchAgglomeration::combineLevels(const label curLevel)
{
    label prevLevel = curLevel - 1;

    // Set the previous level nCells to the current
    nFaces_[prevLevel] = nFaces_[curLevel];

    // Map the restrictAddressing from the coarser level into the previous
    // finer level

    const labelList& curResAddr = restrictAddressing_[curLevel];
    labelList& prevResAddr = restrictAddressing_[prevLevel];

    forAll(prevResAddr, i)
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    // Delete the restrictAddressing for the coarser level
    restrictAddressing_.set(curLevel, nullptr);

    patchLevels_.set(prevLevel, patchLevels_.set(curLevel, nullptr));
}


// ************************************************************************* //
