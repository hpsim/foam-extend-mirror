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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamPointFields_H
#define vtkPVFoamPointFields_H

// Foam includes
#include "interpolatePointToCell.H"

#include "vtkFoamTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPVFoam::convertPointFields
(
    const fvMesh& mesh,
    const pointMesh& pMesh,
    const IOobjectList& objects,
    vtkMultiBlockDataSet* output
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAllConstIter(IOobjectList, objects, iter)
    {
        const word& fieldName = iter()->name();
        // restrict to this GeometricField<Type, ...>
        if
        (
            iter()->headerClassName()
         != GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
        {
            continue;
        }

        if (debug)
        {
            Info<< "Foam::vtkPVFoam::convertPointFields : "
                << fieldName << endl;
        }

        GeometricField<Type, pointPatchField, pointMesh> ptf
        (
            *iter(),
            pMesh
        );


        // Convert activated internalMesh regions
        convertPointFieldBlock
        (
            ptf,
            output,
            partInfoVolume_,
            regionPolyDecomp_
        );

        // Convert activated cellZones
        convertPointFieldBlock
        (
            ptf,
            output,
            partInfoCellZones_,
            zonePolyDecomp_
        );

        // Convert activated cellSets
        convertPointFieldBlock
        (
            ptf,
            output,
            partInfoCellSets_,
            csetPolyDecomp_
        );


        //
        // Convert patches - if activated
        //
        for
        (
            int partId = partInfoPatches_.start();
            partId < partInfoPatches_.end();
            ++partId
        )
        {
            const word  patchName = getPartName(partId);
            const label datasetNo = partDataset_[partId];
            const label   patchId = patches.findPatchID(patchName);

            if (!partStatus_[partId] || datasetNo < 0 || patchId < 0)
            {
                continue;
            }

            convertPatchPointField
            (
                fieldName,
                ptf.boundaryField()[patchId].patchInternalField()(),
                output,
                partInfoPatches_,
                datasetNo
            );
        }
    }
}


template<class Type>
void Foam::vtkPVFoam::convertPointFieldBlock
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const List<polyDecomp>& decompLst
)
{
   for (int partId = selector.start(); partId < selector.end(); ++partId)
   {
       const label datasetNo = partDataset_[partId];

       if (datasetNo >= 0 && partStatus_[partId])
       {
           convertPointField
           (
               ptf,
               GeometricField<Type, fvPatchField, volMesh>::null(),
               output,
               selector,
               datasetNo,
               decompLst[datasetNo]
           );
       }
   }
}


template<class Type>
void Foam::vtkPVFoam::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const label datasetNo,
    const polyDecomp& decomp
)
{
    const label nComp = pTraits<Type>::nComponents;
    const labelList& addPointCellLabels = decomp.addPointCellLabels();
    const labelList& pointMap = decomp.pointMap();

    // use a pointMap or address directly into mesh
    label nPoints;
    if (pointMap.size())
    {
        nPoints = pointMap.size();
    }
    else
    {
        nPoints = ptf.size();
    }

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(nPoints + addPointCellLabels.size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*(nPoints + addPointCellLabels.size()));
    pointData->SetName(ptf.name().c_str());

    if (debug)
    {
        Info<< "convert convertPointField: "
            << ptf.name()
            << " size = " << nPoints
            << " nComp=" << nComp
            << " nTuples = " << (nPoints + addPointCellLabels.size())
            <<  endl;
    }

    float vec[nComp];

    if (pointMap.size())
    {
        forAll(pointMap, i)
        {
            const Type& t = ptf[pointMap[i]];
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
            vtkFoamTupleRemap<Type>(vec);

            pointData->InsertTuple(i, vec);
        }
    }
    else
    {
        forAll(ptf, i)
        {
            const Type& t = ptf[i];
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
            vtkFoamTupleRemap<Type>(vec);

            pointData->InsertTuple(i, vec);
        }
    }

    // continue insertion from here
    label i = nPoints;

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels, apI)
        {
            const Type& t = tf[addPointCellLabels[apI]];
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
            vtkFoamTupleRemap<Type>(vec);

            pointData->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels, apI)
        {
            Type t = interpolatePointToCell(ptf, addPointCellLabels[apI]);
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
            vtkFoamTupleRemap<Type>(vec);

            pointData->InsertTuple(i++, vec);
        }
    }

    vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    )   ->GetPointData()
        ->AddArray(pointData);

    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
