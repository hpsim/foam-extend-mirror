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
#include "faceSet.H"
#include "pointSet.H"
#include "vtkPVFoamPoints.H"

// VTK includes
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::faceSetVTKMesh
(
    const fvMesh& mesh,
    const faceSet& fSet
)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::faceSetVTKMesh" << endl;
        printMemory();
    }

    // Construct primitivePatch of faces in fSet.

    const faceList& meshFaces = mesh.faces();
    faceList patchFaces(fSet.size());
    label faceI = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        patchFaces[faceI++] = meshFaces[iter.key()];
    }
    primitiveFacePatch p(patchFaces, mesh.points());


    // The balance of this routine should be identical to patchVTKMesh

    // Convert Foam mesh vertices to VTK
    const pointField& points = p.localPoints();

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
        Info<< "<end> Foam::vtkPVFoam::faceSetVTKMesh" << endl;
        printMemory();
    }

    return vtkmesh;
}


vtkPolyData* Foam::vtkPVFoam::pointSetVTKMesh
(
    const fvMesh& mesh,
    const pointSet& pSet
)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::pointSetVTKMesh" << endl;
        printMemory();
    }

    const pointField& meshPoints = mesh.points();

    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->Allocate( pSet.size() );

    forAllConstIter(pointSet, pSet, iter)
    {
        vtkPVFoamInsertNextPoint(vtkpoints, meshPoints[iter.key()]);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::pointSetVTKMesh" << endl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
