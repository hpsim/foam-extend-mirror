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

#include "meshToMesh.H"
#include "processorFvPatch.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(meshToMesh, 0);

const Foam::debug::tolerancesSwitch
meshToMesh::directHitTol
(
    "meshToMeshdirectHitTol",
    1e-5
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshToMesh::meshToMesh
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatchNames
)
:
    fromMesh_(meshFrom),
    toMesh_(meshTo),
    patchMap_(patchMap),
    cellAddressing_(toMesh_.nCells()),
    boundaryAddressing_(toMesh_.boundaryMesh().size()),
    inverseDistanceWeightsPtr_(nullptr)
{
    forAll (fromMesh_.boundaryMesh(), patchi)
    {
        fromMeshPatches_.insert
        (
            fromMesh_.boundaryMesh()[patchi].name(),
            patchi
        );
    }

    forAll (toMesh_.boundaryMesh(), patchi)
    {
        toMeshPatches_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            patchi
        );
    }

    forAll (cuttingPatchNames, i)
    {
        if (toMeshPatches_.found(cuttingPatchNames[i]))
        {
            cuttingPatches_.insert
            (
                cuttingPatchNames[i],
                toMeshPatches_.find(cuttingPatchNames[i])()
            );
        }
        else
        {
            WarningIn
            (
                "meshToMesh::meshToMesh"
                "(const fvMesh& meshFrom, const fvMesh& meshTo,"
                "const HashTable<word>& patchMap,"
                "const wordList& cuttingPatchNames)"
            )   << "Cannot find cutting-patch " << cuttingPatchNames[i]
                << " in destination mesh" << endl;
        }
    }

    forAll (toMesh_.boundaryMesh(), patchi)
    {
        // Add the processor patches in the toMesh to the cuttingPatches list
        if (isA<processorPolyPatch>(toMesh_.boundaryMesh()[patchi]))
        {
            cuttingPatches_.insert
            (
                toMesh_.boundaryMesh()[patchi].name(),
                patchi
            );
        }
    }

    calcAddressing();
}


meshToMesh::meshToMesh
(
    const fvMesh& meshFrom,
    const fvMesh& meshTo
)
:
    fromMesh_(meshFrom),
    toMesh_(meshTo),
    cellAddressing_(toMesh_.nCells()),
    boundaryAddressing_(toMesh_.boundaryMesh().size()),
    inverseDistanceWeightsPtr_(nullptr)
{
    // check whether both meshes have got the same number
    // of boundary patches
    if (fromMesh_.boundary().size() != toMesh_.boundary().size())
    {
        FatalErrorIn
        (
            "meshToMesh::meshToMesh"
            "(const fvMesh& meshFrom, const fvMesh& meshTo)"
        )   << "Incompatible meshes: different number of patches, "
            << "fromMesh = " << fromMesh_.boundary().size()
            << ", toMesh = " << toMesh_.boundary().size()
            << exit(FatalError);
    }

    forAll (fromMesh_.boundaryMesh(), patchi)
    {
        if
        (
            fromMesh_.boundaryMesh()[patchi].name()
         != toMesh_.boundaryMesh()[patchi].name()
        )
        {
            FatalErrorIn
            (
                "meshToMesh::meshToMesh"
                "(const fvMesh& meshFrom, const fvMesh& meshTo)"
            )   << "Incompatible meshes: different patch names for patch "
                << patchi
                << ", fromMesh = " << fromMesh_.boundary()[patchi].name()
                << ", toMesh = " << toMesh_.boundary()[patchi].name()
                << exit(FatalError);
        }

        if
        (
            fromMesh_.boundaryMesh()[patchi].type()
         != toMesh_.boundaryMesh()[patchi].type()
        )
        {
            FatalErrorIn
            (
                "meshToMesh::meshToMesh"
                "(const fvMesh& meshFrom, const fvMesh& meshTo)"
            )   << "Incompatible meshes: different patch types for patch "
                << patchi
                << ", fromMesh = " << fromMesh_.boundary()[patchi].type()
                << ", toMesh = " << toMesh_.boundary()[patchi].type()
                << exit(FatalError);
        }

        fromMeshPatches_.insert
        (
            fromMesh_.boundaryMesh()[patchi].name(),
            patchi
        );

        toMeshPatches_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            patchi
        );

        patchMap_.insert
        (
            toMesh_.boundaryMesh()[patchi].name(),
            fromMesh_.boundaryMesh()[patchi].name()
        );
    }

    calcAddressing();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshToMesh::~meshToMesh()
{
    deleteDemandDrivenData(inverseDistanceWeightsPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
