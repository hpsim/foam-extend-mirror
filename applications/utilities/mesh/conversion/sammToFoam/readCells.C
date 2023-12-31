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
    Create intermediate mesh from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void sammMesh::addRegularCell
(
    const labelList& labels,
    const label nCreatedCells
)
{
    // Momory management
    static labelList labelsHex(8);
    static labelList labelsWedge(7);
    static labelList labelsPrism(6);
    static labelList labelsPyramid(5);
    static labelList labelsTet(4);
    static labelList labelsTetWedge(5);

    if      // Tetrahedron
    (
        labels[2] == labels[3]
     && labels[4] == labels[5]
     && labels[5] == labels[6]
     && labels[6] == labels[7]
    )
    {
        labelsTet[0] = labels[0];
        labelsTet[1] = labels[1];
        labelsTet[2] = labels[2];
        labelsTet[3] = labels[4];
        cellShapes_[nCreatedCells] = cellShape(*tetPtr_, labelsTet);
    }

    else if // Square-based pyramid
    (
        labels[4] == labels[5]
     && labels[5] == labels[6]
     && labels[6] == labels[7]
    )
    {
        labelsPyramid[0] = labels[0];
        labelsPyramid[1] = labels[1];
        labelsPyramid[2] = labels[2];
        labelsPyramid[3] = labels[3];
        labelsPyramid[4] = labels[4];
        cellShapes_[nCreatedCells] = cellShape(*pyrPtr_, labelsPyramid);
    }

    else if // Tet Wedge
    (
        labels[2] == labels[3]
     && labels[4] == labels[5]
     && labels[6] == labels[7]
    )
    {
        labelsTetWedge[0] = labels[0];
        labelsTetWedge[1] = labels[1];
        labelsTetWedge[2] = labels[2];
        labelsTetWedge[3] = labels[4];
        labelsTetWedge[4] = labels[6];
        cellShapes_[nCreatedCells] = cellShape(*tetWedgePtr_, labelsTetWedge);
    }

    else if // Triangular prism
    (
        labels[2] == labels[3]
     && labels[6] == labels[7]
    )
    {
        labelsPrism[0] = labels[0];
        labelsPrism[1] = labels[1];
        labelsPrism[2] = labels[2];
        labelsPrism[3] = labels[4];
        labelsPrism[4] = labels[5];
        labelsPrism[5] = labels[6];
        cellShapes_[nCreatedCells] = cellShape(*prismPtr_, labelsPrism);
    }

    else if // Wedge
    (
        labels[4] == labels[7]
    )
    {
        labelsWedge[0] = labels[7];
        labelsWedge[1] = labels[6];
        labelsWedge[2] = labels[5];
        labelsWedge[3] = labels[3];
        labelsWedge[4] = labels[2];
        labelsWedge[5] = labels[1];
        labelsWedge[6] = labels[0];
        cellShapes_[nCreatedCells] = cellShape(*wedgePtr_, labelsWedge);
    }

    else    // Hex
    {
        labelsHex[0] = labels[0];
        labelsHex[1] = labels[1];
        labelsHex[2] = labels[2];
        labelsHex[3] = labels[3];
        labelsHex[4] = labels[4];
        labelsHex[5] = labels[5];
        labelsHex[6] = labels[6];
        labelsHex[7] = labels[7];
        cellShapes_[nCreatedCells] = cellShape(*hexPtr_, labelsHex);
    }
}


void sammMesh::addSAMMcell
(
    const label typeFlag,
    const labelList& globalLabels,
    const label nCreatedCells
)
{

    // grab the shape from the table
    if (!sammShapeLookup[typeFlag] || !sammAddressingTable[typeFlag])
    {
        FatalErrorIn
        (
            "sammMesh::addRegularCell(const labelList& labels, "
            "const label nCreatedCells)"
        )   << "SAMM type " << typeFlag << " has no registered label. BUG!"
            << abort(FatalError);
    }

     const cellModel& curModel = *(sammShapeLookup[typeFlag]);

    // get reference to the addressing list
    const label* addressing = sammAddressingTable[typeFlag];

    // make a list of labels
    labelList sammCellLabels(curModel.nPoints(), -1);

    forAll (sammCellLabels, labelI)
    {
        sammCellLabels[labelI] = globalLabels[addressing[labelI]];
    }

    cellShapes_[nCreatedCells] = cellShape(curModel, sammCellLabels);
}


void sammMesh::readCells()
{
    label nCells = 0;
    label maxLabel = -1;

    fileName cellsFileName(casePrefix_ + ".cel");

    {
        IFstream cellsFile(cellsFileName);

        if (cellsFile.good())
        {
            label lineLabel, cellLabel = -1, pointLabel, regionLabel, typeFlag;

            maxLabel = -1;
            while (!(cellsFile >> lineLabel).eof())
            {
                maxLabel = max(maxLabel, lineLabel);
                for (int i=0; i<8; i++)
                {
                    cellsFile >> pointLabel;
                }

                cellsFile >> regionLabel;
                cellsFile >> typeFlag;

                if (lineLabel != cellLabel)
                {
                    cellLabel = lineLabel;
                    nCells++;
                }
            }
        }
        else
        {
            FatalErrorIn("sammMesh::readCells()")
                << "Cannot read file "
                << cellsFileName
                << abort(FatalError);
        }
    }

    Info<< "Number of cells = " << nCells << endl << endl;

    cellShapes_.setSize(nCells);

    starCellLabelLookup_.setSize(maxLabel+1);

    // reset point labels to invalid value
    forAll (starCellLabelLookup_, i)
    {
        starCellLabelLookup_[i] = -1;
    }


    if (nCells > 0)
    {
        IFstream cellsFile(cellsFileName);

        labelList labels(24, label(-1));
        label lineLabel, sammLabel, regionLabel, typeFlag;

        for (label cellI = 0; cellI < nCells; cellI++)
        {
            label nLabels = 0;

            bool addOnToCell = false;

            do
            {
                if (nLabels > 24)
                {
                    FatalErrorIn("sammMesh::readCells()")
                        << "Unknown SAMM cell. "
                        << "More than 24 vertices"
                        << abort(FatalError);
                }

                if ((cellsFile >> lineLabel).eof())
                {
                    FatalErrorIn("sammMesh::readCells()")
                        << "Reached end of cells file before "
                        << "all cells are read in."
                        << abort(FatalError);
                }

                // prepare for possible continuation
                nLabels += 8;

                for (int i=nLabels-8; i<nLabels; i++)
                {
                    cellsFile >> sammLabel;

                    if (sammLabel != 0)
                    {
                        // Convert Samm vertex number to point label
                        labels[i] = starPointLabelLookup_[sammLabel];

                        if (labels[i] < 0)
                        {
                            Info<< "Cell file not consistent with vertex file. "
                                << "Samm vertex number " << sammLabel
                                << " does not exist\n";
                        }
                    }
                    else
                    {
                        labels[i] = -1;
                    }
                }

                cellsFile >> regionLabel;
                cellsFile >> typeFlag;

                // check for continuation line
                if (!addOnToCell && typeFlag == 255)
                {
                    addOnToCell = true;
                }
                else
                {
                    addOnToCell = false;
                }

            } while (typeFlag == -1 || addOnToCell);

            starCellLabelLookup_[lineLabel] = cellI;

            if (nLabels == 8)
            {
                addRegularCell(labels, cellI);
            }
            else
            {
                addSAMMcell(typeFlag, labels, cellI);
            }
        }
    }
    else
    {
        FatalErrorIn("sammMesh::readCells()")
            << "No cells in file "
            << cellsFileName
            << abort(FatalError);
    }
}


// ************************************************************************* //
