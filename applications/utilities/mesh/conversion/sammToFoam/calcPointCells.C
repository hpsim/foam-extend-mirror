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
    Create intermediate mesh files from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sammMesh::calcPointCells() const
{
    const static label UNIT_POINT_CELLS = 12;

    if (pointCellsPtr_)
    {
        FatalErrorIn("sammMesh::calcPointCells()")
            << "PointCells already calculated"
            << abort(FatalError);
    }


    pointCellsPtr_ = new labelListList(points_.size());

    labelListList& pc = *pointCellsPtr_;

    forAll(pc, i)
    {
        pc[i].setSize(UNIT_POINT_CELLS);
    }

    // Initialise the list of labels which will hold the count the
    // actual number of cells per point during the analysis
    labelList cellCount(points_.size());

    forAll(cellCount, i)
    {
        cellCount[i] = 0;
    }

    // Note. Unlike the standard point-cell algorithm, which asks the cell for
    // the supporting point labels, we need to work based on the cell faces.
    // This is because some of the faces for meshes with arbitrary interfaces
    // do not come from the cell shape, but from the slaves of the coupled
    // match. It is also adventageous to remove the duplicates from the
    // point-cell addressing, because this removes a lot of waste later.
    //

    // For each cell
    forAll(cellShapes_, cellI)
    {
        const faceList& faces = cellFaces_[cellI];

        forAll (faces, i)
        {
            // For each vertex
            const labelList& labels = faces[i];

            forAll(labels, j)
            {
                // Set working point label
                label curPoint = labels[j];
                labelList& curPointCells = pc[curPoint];
                label curCount = cellCount[curPoint];

                // check if the cell has been added before
                bool found = false;

                for (label f = 0; f < curCount; f++)
                {
                    if (curPointCells[f] == cellI)
                    {
                        found = true;

                        break;
                    }
                }

                if (!found)
                {

                    // If the list of pointCells is not big enough, double it
                    if (curPointCells.size() <= curCount)
                    {
                        curPointCells.setSize(curPointCells.size()*2);
                    }

                    // Enter the cell label in the point's cell list
                    curPointCells[curCount] = cellI;

                    // Increment the cell count for the point addressed
                    cellCount[curPoint]++;
                }
            }
        }
    }

    // Finally, truncate the lists made to their active size
    forAll(pc, i)
    {
        pc[i].setSize(cellCount[i]);
    }
}


const labelListList& sammMesh::pointCells() const
{
    if (!pointCellsPtr_)
    {
        calcPointCells();
    }

    return *pointCellsPtr_;
}


// ************************************************************************* //
