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

#include "writePatch.H"
#include "OFstream.H"
#include "writeFuns.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writePatch
(
    const bool binary,
    const word& setName,
    const primitiveFacePatch& fp,
    const word& fieldName,
    labelList& fieldValues,
    const fileName& fileName
)
{
    std::ofstream pStream(fileName.c_str());

    pStream
        << "# vtk DataFile Version 2.0" << std::endl
        << setName << std::endl;
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


    // Write points and faces as polygons

    pStream << "POINTS " << fp.nPoints() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*fp.nPoints());

    writeFuns::insert(fp.localPoints(), ptField);

    writeFuns::write(pStream, binary, ptField);


    label nFaceVerts = 0;

    forAll(fp.localFaces(), faceI)
    {
        nFaceVerts += fp.localFaces()[faceI].size() + 1;
    }
    pStream << "POLYGONS " << fp.size() << ' ' << nFaceVerts
        << std::endl;


    dynamicLabelList vertLabels(nFaceVerts);

    forAll(fp.localFaces(), faceI)
    {
        const face& f = fp.localFaces()[faceI];

        vertLabels.append(f.size());

        writeFuns::insert(f, vertLabels);
    }
    writeFuns::write(pStream, binary, vertLabels);


    //-----------------------------------------------------------------
    //
    // Write data
    //
    //-----------------------------------------------------------------

    // Write faceID

    pStream
        << "CELL_DATA " << fp.size() << std::endl
        << "FIELD attributes 1" << std::endl;

    // Cell ids first
    pStream << fieldName << " 1 " << fp.size() << " int" << std::endl;

    writeFuns::write(pStream, binary, fieldValues);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
