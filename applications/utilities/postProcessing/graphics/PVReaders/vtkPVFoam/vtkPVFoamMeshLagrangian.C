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
#include "CloudTemplate.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "passiveParticle.H"
#include "vtkPVFoamPoints.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::lagrangianVTKMesh
(
    const fvMesh& mesh,
    const word& cloudName
)
{
    vtkPolyData* vtkmesh = nullptr;

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::lagrangianVTKMesh - timePath "
            << mesh.time().timePath()/cloud::prefix/cloudName << endl;
        printMemory();
    }


    // the region name is already in the mesh db
    IOobjectList sprayObjs
    (
        mesh,
        mesh.time().timeName(),
        cloud::prefix/cloudName
    );

    IOobject* positionsPtr = sprayObjs.lookup("positions");
    if (positionsPtr)
    {
        Cloud<passiveParticle> parcels(mesh, cloudName, false);

        if (debug)
        {
            Info<< "cloud with " << parcels.size() << " parcels" << endl;
        }

        vtkmesh = vtkPolyData::New();
        vtkPoints* vtkpoints = vtkPoints::New();
        vtkCellArray* vtkcells = vtkCellArray::New();

        vtkpoints->Allocate( parcels.size() );
        vtkcells->Allocate( parcels.size() );

        vtkIdType particleId = 0;
        forAllConstIter(Cloud<passiveParticle>, parcels, iter)
        {
            vtkPVFoamInsertNextPoint(vtkpoints, iter().position());

            vtkcells->InsertNextCell(1, &particleId);
            particleId++;
        }

        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();

        vtkmesh->SetVerts(vtkcells);
        vtkcells->Delete();
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::lagrangianVTKMesh" << endl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
