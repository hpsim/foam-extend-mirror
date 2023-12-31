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

#include "pointFieldReconstructor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointFieldReconstructor::pointFieldReconstructor
(
    const pointMesh& mesh,
    const PtrList<pointMesh>& procMeshes,
    const PtrList<labelIOList>& pointProcAddressing,
    const PtrList<labelIOList>& boundaryProcAddressing
)
:
    mesh_(mesh),
    procMeshes_(procMeshes),
    pointProcAddressing_(pointProcAddressing),
    boundaryProcAddressing_(boundaryProcAddressing),
    patchPointAddressing_(procMeshes.size())
{
    // Inverse-addressing of the patch point labels.
    labelList pointMap(mesh_.size(), -1);

    // Create the pointPatch addressing
    forAll(procMeshes_, proci)
    {
        const pointMesh& procMesh = procMeshes_[proci];

        patchPointAddressing_[proci].setSize(procMesh.boundary().size());

        forAll(procMesh.boundary(), patchi)
        {
            if (boundaryProcAddressing_[proci][patchi] >= 0)
            {
                labelList& procPatchAddr = patchPointAddressing_[proci][patchi];
                procPatchAddr.setSize(procMesh.boundary()[patchi].size(), -1);

                const labelList& patchPointLabels =
                    mesh_.boundary()[boundaryProcAddressing_[proci][patchi]]
                    .meshPoints();

                // Create the inverse-addressing of the patch point labels.
                forAll (patchPointLabels, pointi)
                {
                    pointMap[patchPointLabels[pointi]] = pointi;
                }

                const labelList& procPatchPoints =
                    procMesh.boundary()[patchi].meshPoints();

                forAll (procPatchPoints, pointi)
                {
                    procPatchAddr[pointi] =
                        pointMap
                        [
                            pointProcAddressing_[proci][procPatchPoints[pointi]]
                        ];
                }

                if (procPatchAddr.size() && min(procPatchAddr) < 0)
                {
                    FatalErrorIn
                    (
                        "pointFieldReconstructor::pointFieldReconstructor"
                        "(\n"
                        "    const pointMesh& mesh,\n"
                        "    const PtrList<pointMesh>& procMeshes,\n"
                        "    const PtrList<labelIOList>& pointProcAddressing,\n"
                        "    const PtrList<labelIOList>& "
                        "boundaryProcAddressing\n"
                        ")"
                    )   << "Incomplete patch point addressing"
                        << abort(FatalError);
                }
            }
        }
    }
}


// ************************************************************************* //
