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

#include "writePointSet.H"
#include "OFstream.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writePointSet
(
    const bool binary,
    const vtkMesh& vMesh,
    const pointSet& set,
    const fileName& fileName
)
{
    std::ofstream pStream(fileName.c_str());

    pStream
        << "# vtk DataFile Version 2.0" << std::endl
        << set.name() << std::endl;
    if (binary)
    {
        pStream << "BINARY" << std::endl;
    }
    else
    {
        pStream << "ASCII" << std::endl;
    }
    pStream << "DATASET POLYDATA" << std::endl;


    //------------------------------------------------------------------
    //
    // Write topology
    //
    //------------------------------------------------------------------


    // Write points

    pStream << "POINTS " << set.size() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*set.size());

    writeFuns::insert
    (
        UIndirectList<point>(vMesh.mesh().points(), set.toc())(),
        ptField
    );

    writeFuns::write(pStream, binary, ptField);


    //-----------------------------------------------------------------
    //
    // Write data
    //
    //-----------------------------------------------------------------

    // Write faceID

    pStream
        << "POINT_DATA " << set.size() << std::endl
        << "FIELD attributes 1" << std::endl;

    // Cell ids first
    pStream << "pointID 1 " << set.size() << " int" << std::endl;

    labelList pointIDs(set.toc());

    writeFuns::write(pStream, binary, pointIDs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
