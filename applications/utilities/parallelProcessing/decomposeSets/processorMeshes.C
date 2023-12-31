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

#include "processorMeshes.H"
#include "foamTime.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorMeshes::readMeshes()
{
    forAll (databases_, procI)
    {
        meshes_.set
        (
            procI,
            new fvMesh
            (
                IOobject
                (
                    meshName_,
                    databases_[procI].timeName(),
                    databases_[procI]
                )
            )
        );
    }
}


void Foam::processorMeshes::readAddressing()
{
    forAll (databases_, procI)
    {
        pointProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "pointProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        faceProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        cellProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        boundaryProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorMeshes::processorMeshes
(
    PtrList<Time>& databases,
    const word& meshName
)
:
    databases_(databases),
    meshName_(meshName),
    meshes_(databases.size()),
    pointProcAddressing_(databases.size()),
    faceProcAddressing_(databases.size()),
    cellProcAddressing_(databases.size()),
    boundaryProcAddressing_(databases.size())
{
    readMeshes();
    readAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState Foam::processorMeshes::readUpdate()
{
    polyMesh::readUpdateState stat = polyMesh::UNCHANGED;

    forAll (databases_, procI)
    {
        // Check if any new meshes need to be read.
        polyMesh::readUpdateState procStat = meshes_[procI].readUpdate();

        // Force changed mesh if the processor mesh holds points. This is
        // essential for changed region meshes and if the current time is
        // equal to the initial time during construction.
        // Pascal Beckstein, 14/Jul/2016
        fileName timePath = meshes_[procI].time().timePath();
        if (isFile(timePath/meshes_[procI].meshDir()/"points"))
        {
            procStat = polyMesh::POINTS_MOVED;
        }

        // Combine into overall mesh change status
        if (stat == polyMesh::UNCHANGED)
        {
            stat = procStat;
        }
        else
        {
            if (stat != procStat)
            {
                FatalErrorIn("processorMeshes::readUpdate()")
                    << "Processor " << procI
                    << " has a different polyMesh at time "
                    << databases_[procI].timeName()
                    << " compared to any previous processors." << nl
                    << "Please check time " << databases_[procI].timeName()
                    << " directories on all processors for consistent"
                    << " mesh files."
                    << exit(FatalError);
            }
        }
    }

    // Reading of meshes removed: readUpdate will do this
    if
    (
        stat == polyMesh::TOPO_CHANGE
     || stat == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        // Reread addressing; meshes are already updated with readUpdate.
        readAddressing();
    }

    return stat;
}


void Foam::processorMeshes::reconstructPoints(fvMesh& mesh)
{
    // We may not use
    //   mesh.movePoints(newPoints);
    //   mesh.write();
    // here, as this will fail for mesh regions containing
    //   directMappedPolyPatches or
    //   directMappedWallPolyPatches
    // and probably other region-coupled patches. E.g. for both above
    // mentioned directMapped patch types we have calls to
    //   directMappedPatchBase::map()
    // in the corresponding initMovePoints(const pointField& p) member
    // function. This however involves parallel communication and it would
    // require all other connected regions to be already reconstructed and
    // registered during reconstruction of the current region!
    // Pascal Beckstein, 14/Jul/2016

    // Create the new points
    pointIOField newPoints
    (
        IOobject
        (
            "points",
            mesh.time().timeName(),
            mesh.meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh.allPoints()
    );

    forAll (meshes_, procI)
    {
        const vectorField& procPoints = meshes_[procI].allPoints();

        // Set the cell values in the reconstructed field

        const labelList& pointProcAddressingI = pointProcAddressing_[procI];

        if (pointProcAddressingI.size() != procPoints.size())
        {
            FatalErrorIn("processorMeshes")
                << "problem :"
                << " pointProcAddressingI:" << pointProcAddressingI.size()
                << " procPoints:" << procPoints.size()
                << abort(FatalError);
        }

        // Only live points carry reconstruction data.  Reconsider
        // HJ, 6/Sep/2009
        for (label pointI = 0; pointI < meshes_[procI].nPoints(); pointI++)
//         forAll(pointProcAddressingI, pointI)
        {
            newPoints[pointProcAddressingI[pointI]] = procPoints[pointI];
        }
    }

    newPoints.write();
    mesh.readUpdate();
}


// ************************************************************************* //
