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
    Create intermediate mesh from PROSTAR files

\*---------------------------------------------------------------------------*/

#include "starMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void starMesh::readPoints(const scalar scaleFactor)
{
    label nPoints = 0;
    label maxLabel = -1;

    fileName pointsFileName(casePrefix_ + ".vrt");

    {
        IFstream pointsFile(pointsFileName);

        // Pass 1: get # points and maximum vertex label

        if (pointsFile.good())
        {
            label pointLabel;
            scalar x, y, z;

            maxLabel = -1;
            while ((pointsFile >> pointLabel).good())
            {
                nPoints++;
                maxLabel = max(maxLabel, pointLabel);
                pointsFile >> x >> y >> z;
            }
        }
        else
        {
            FatalErrorIn("starMesh::readPoints()")
                << "Cannot read file " << pointsFileName
                << abort(FatalError);
        }
    }

    Info<< "Number of points = " << nPoints << endl << endl;

    points_.setSize(nPoints);
    starPointID_.setSize(nPoints);

    // Reset STAR point ID, just in case
    starPointID_ = -1;

    starPointLabelLookup_.setSize(maxLabel+1);

    // reset point labels to invalid value
    starPointLabelLookup_ = -1;

    if (nPoints > 0)
    {
        // Pass 2: construct pointlist and conversion table
        // from Star vertex numbers to Foam pointLabels

        IFstream pointsFile(pointsFileName);
        label pointLabel;

        forAll(points_, p)
        {
            pointsFile
                >> pointLabel
                >> points_[p].x()
                >> points_[p].y()
                >> points_[p].z();

            starPointID_[p] = pointLabel;
            starPointLabelLookup_[pointLabel] = p;
        }

        if (scaleFactor > 1.0 + SMALL || scaleFactor < 1.0 - SMALL)
        {
            points_ *= scaleFactor;
        }
    }
    else
    {
        FatalError
            << "void starMesh::readPoints() : "
            << "no points in file "
            << pointsFileName
            << abort(FatalError);
    }
}


// ************************************************************************* //
