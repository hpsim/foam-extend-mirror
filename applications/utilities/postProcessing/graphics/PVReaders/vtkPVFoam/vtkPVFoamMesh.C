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

#include "vtkPVFoam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fvMeshSubset.H"
#include "vtkPVFoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertMeshVolume
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoVolume_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // resize for decomposed polyhedra
    regionPolyDecomp_.setSize(selector.size());

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshVolume" << endl;
        printMemory();
    }

    // Convert the internalMesh
    // this looks like more than one part, but it isn't
    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word partName = "internalMesh";

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            mesh,
            regionPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, selector, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshVolume" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoLagrangian_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshLagrangian" << endl;
        printMemory();
    }

    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word cloudName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkPolyData* vtkmesh = lagrangianVTKMesh(mesh, cloudName);

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, selector, datasetNo, cloudName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshLagrangian" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPatches
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoPatches_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshPatches" << endl;
        printMemory();
    }

    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word patchName = getPartName(partId);
        const label  patchId = patches.findPatchID(patchName);

        if (!partStatus_[partId] || patchId < 0)
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for patch[" << patchId <<"] "
                << patchName  << endl;
        }

        vtkPolyData* vtkmesh = patchVTKMesh(patches[patchId]);

        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, selector, datasetNo, patchName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshPatches" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoCellZones_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // resize for decomposed polyhedra
    zonePolyDecomp_.setSize(selector.size());

    if (!selector.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshCellZones" << endl;
        printMemory();
    }

    const cellZoneMesh& zMesh = mesh.cellZones();
    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (!partStatus_[partId] || zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for cellZone[" << zoneId << "] "
                << zoneName << endl;
        }

        fvMeshSubset subsetMesh
        (
            IOobject
            (
                "set",
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        subsetMesh.setLargeCellSubset(labelHashSet(zMesh[zoneId]));

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetMesh.subMesh(),
            zonePolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetMesh.cellMap(),
                zonePolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetMesh.cellMap(),
                zonePolyDecomp_[datasetNo].addPointCellLabels()
            );

            // copy pointMap as well, otherwise pointFields fail
            zonePolyDecomp_[datasetNo].pointMap() = subsetMesh.pointMap();

            AddToBlock(output, vtkmesh, selector, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshCellZones" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoCellSets_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // resize for decomposed polyhedra
    csetPolyDecomp_.setSize(selector.size());

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshCellSets" << endl;
        printMemory();
    }

    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for cellSet=" << partName << endl;
        }

        const cellSet cSet(mesh, partName);

        fvMeshSubset subsetMesh
        (
            IOobject
            (
                "set",
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        subsetMesh.setLargeCellSubset(cSet);

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetMesh.subMesh(),
            csetPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetMesh.cellMap(),
                csetPolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetMesh.cellMap(),
                csetPolyDecomp_[datasetNo].addPointCellLabels()
            );

            // copy pointMap as well, otherwise pointFields fail
            csetPolyDecomp_[datasetNo].pointMap() = subsetMesh.pointMap();

            AddToBlock(output, vtkmesh, selector, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshCellSets" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoFaceZones_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (!selector.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshFaceZones" << endl;
        printMemory();
    }

    const faceZoneMesh& zMesh = mesh.faceZones();
    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (!partStatus_[partId] || zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTKmesh for faceZone[" << zoneId << "] "
                << zoneName << endl;
        }

        vtkPolyData* vtkmesh = faceZoneVTKMesh(mesh, zMesh[zoneId]);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, selector, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshFaceZones" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoFaceSets_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshFaceSets" << endl;
        printMemory();
    }

    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for faceSet=" << partName << endl;
        }

        const faceSet fSet(mesh, partName);

        vtkPolyData* vtkmesh = faceSetVTKMesh(mesh, fSet);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, selector, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshFaceSets" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoPointZones_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshPointZones" << endl;
        printMemory();
    }

    if (selector.size())
    {
        const pointZoneMesh& zMesh = mesh.pointZones();
        for (int partId = selector.start(); partId < selector.end(); ++partId)
        {
            word zoneName = getPartName(partId);
            label zoneId = zMesh.findZoneID(zoneName);

            if (!partStatus_[partId] || zoneId < 0)
            {
                continue;
            }

            vtkPolyData* vtkmesh = pointZoneVTKMesh(mesh, zMesh[zoneId]);
            if (vtkmesh)
            {
                AddToBlock(output, vtkmesh, selector, datasetNo, zoneName);
                vtkmesh->Delete();

                partDataset_[partId] = datasetNo++;
            }
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshPointZones" << endl;
        printMemory();
    }
}



void Foam::vtkPVFoam::convertMeshPointSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoPointSets_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::convertMeshPointSets" << endl;
        printMemory();
    }

    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for pointSet=" << partName << endl;
        }

        const pointSet pSet(mesh, partName);

        vtkPolyData* vtkmesh = pointSetVTKMesh(mesh, pSet);
        if (vtkmesh)
        {
            AddToBlock(output, vtkmesh, selector, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::convertMeshPointSets" << endl;
        printMemory();
    }
}

// ************************************************************************* //
