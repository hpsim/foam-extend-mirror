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

#include "reconstructLagrangian.H"
#include "labelIOList.H"
#include "CloudTemplate.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::reconstructLagrangianPositions
(
    const polyMesh& mesh,
    const word& cloudName,
    const PtrList<fvMesh>& meshes,
    const PtrList<labelIOList>& faceProcAddressing,
    const PtrList<labelIOList>& cellProcAddressing
)
{
    Cloud<passiveParticle> lagrangianPositions
    (
        mesh,
        cloudName,
        IDLList<passiveParticle>()
    );

    forAll (meshes, i)
    {
        if (meshes.set(i))
        {
            const labelList& cellMap = cellProcAddressing[i];

            Cloud<passiveParticle> lpi(meshes[i], cloudName, false);

            forAllIter (Cloud<passiveParticle>, lpi, iter)
            {
                const passiveParticle& ppi = iter();

                lagrangianPositions.append
                (
                    new passiveParticle
                    (
                        lagrangianPositions,
                        ppi.position(),
                        cellMap[ppi.cell()]
                    )
                );
            }
        }
    }

    IOPosition<passiveParticle>(lagrangianPositions).write();
}


// ************************************************************************* //
