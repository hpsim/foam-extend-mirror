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

#include "writeSurfFields.H"
#include "OFstream.H"
#include "floatScalar.H"
#include "writeFuns.H"
#include "emptyFvsPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writeSurfFields
(
    const bool binary,
    const vtkMesh& vMesh,
    const fileName& fileName,
    const PtrList<surfaceVectorField>& surfVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();

    std::ofstream str(fileName.c_str());

    str << "# vtk DataFile Version 2.0" << std::endl
        << "surfaceFields" << std::endl;

    if (binary)
    {
        str << "BINARY" << std::endl;
    }
    else
    {
        str << "ASCII" << std::endl;
    }
    str << "DATASET POLYDATA" << std::endl;

    const pointField& fc = mesh.faceCentres();

    str << "POINTS " << mesh.nFaces() << " float" << std::endl;

    DynamicList<floatScalar> pField(3*mesh.nFaces());

    for (label faceI = 0; faceI < mesh.nFaces(); faceI++)
    {
        writeFuns::insert(fc[faceI], pField);
    }

    writeFuns::write(str, binary, pField);

    str << "POINT_DATA " << mesh.nFaces() << std::endl
        << "FIELD attributes " << surfVectorFields.size() << std::endl;

    // surfVectorFields
    forAll(surfVectorFields, fieldI)
    {
        const surfaceVectorField& svf = surfVectorFields[fieldI];

        str << svf.name() << " 3 "
            << mesh.nFaces() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*mesh.nFaces());

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            writeFuns::insert(svf[faceI], fField);
        }

        forAll(svf.boundaryField(), patchI)
        {
            const fvsPatchVectorField& pf = svf.boundaryField()[patchI];

            const fvPatch& pp = mesh.boundary()[patchI];

            if (isA<emptyFvsPatchVectorField>(pf))
            {
                // Note: loop over polypatch size, not fvpatch size.
                forAll(pp.patch(), i)
                {
                    writeFuns::insert(vector::zero, fField);
                }
            }
            else
            {
                forAll(pf, i)
                {
                    writeFuns::insert(pf[i], fField);
                }
            }
        }

        writeFuns::write(str, binary, fField);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
