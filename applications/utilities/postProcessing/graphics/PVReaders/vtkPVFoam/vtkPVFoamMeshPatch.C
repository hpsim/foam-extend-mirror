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
#include "polyPatch.H"
#include "primitivePatch.H"
#include "vtkPVFoamPoints.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::patchVTKMesh
(
    const polyPatch& p
)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::patchVTKMesh - " << p.name() << endl;
        printMemory();
    }

    // Convert Foam mesh vertices to VTK
    const Foam::pointField& points = p.localPoints();

    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->Allocate( points.size() );
    forAll(points, i)
    {
        vtkPVFoamInsertNextPoint(vtkpoints, points[i]);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();


    // Add faces as polygons
    const faceList& faces = p.localFaces();

    vtkCellArray* vtkcells = vtkCellArray::New();
    vtkcells->Allocate( faces.size() );
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];
        vtkIdType nodeIds[f.size()];

        forAll(f, fp)
        {
            nodeIds[fp] = f[fp];
        }
        vtkcells->InsertNextCell(f.size(), nodeIds);
    }

    vtkmesh->SetPolys(vtkcells);
    vtkcells->Delete();

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::patchVTKMesh - " << p.name() << endl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
