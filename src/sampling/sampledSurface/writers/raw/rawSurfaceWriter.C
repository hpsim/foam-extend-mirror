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

#include "rawSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"

#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(rawSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::rawSurfaceWriter::writeLocation
(
    Ostream& os,
    const pointField& points,
    const label pointI
)
{
    const point& pt = points[pointI];
    os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ';
}


inline void Foam::rawSurfaceWriter::writeLocation
(
    Ostream& os,
    const pointField& points,
    const faceList& faces,
    const label faceI
)
{
    const point& ct = faces[faceI].centre(points);
    os  << ct.x() << ' ' << ct.y() << ' ' << ct.z() << ' ';
}


namespace Foam
{
    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<scalar>& values
    )
    {
        os  << values.size() << nl
            << "#  x  y  z  " << fieldName << nl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<vector>& values
    )
    {
        os  << values.size() << nl
            << "#  x  y  z  "
            << fieldName << "_x  "
            << fieldName << "_y  "
            << fieldName << "_z  "
            << endl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<sphericalTensor>& values
    )
    {
        os  << values.size() << nl
            << "#  ii  "
            << fieldName << "_ii" << nl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<symmTensor>& values
    )
    {
        os  << values.size() << nl
            << "#  xx  xy  xz  yy  yz ";
        for (int i=0; i<6; ++i)
        {
            os  << fieldName << "_" << i << "  ";
        }
        os  << endl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<symmTensor4thOrder>& values
    )
    {
        os  << values.size() << nl
            << "#  xxxx  xxyy  xxzz  yyyy  yyzz zzzz xyxy yzyzy zxzx";
        for (int i=0; i<9; ++i)
        {
            os  << fieldName << "_" << i << "  ";
        }
        os  << endl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<diagTensor>& values
    )
    {
        os  << values.size() << nl
            << "#  xx  yy  zz ";
        for (int i=0; i<3; ++i)
        {
            os  << fieldName << "_" << i << "  ";
        }
        os  << endl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<tensor>& values
    )
    {
        os  << values.size() << nl
            << "#  xx  xy  xz  yx  yy  yz  zx  zy  zz";
        for (int i=0; i<9; ++i)
        {
            os  << fieldName << "_" << i << "  ";
        }
        os  << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const scalar& v
    )
    {
        os  << v << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const vector& v
    )
    {
        os  << v[0] << ' ' << v[1] << ' ' << v[2] << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const sphericalTensor& v
    )
    {
        os  << v[0] << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const symmTensor& v
    )
    {
        os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
            << v[3] << ' ' << v[4] << ' ' << v[5] << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const symmTensor4thOrder& v
    )
    {
        os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
            << v[3] << ' ' << v[4] << ' ' << v[5] << ' '
            << v[6] << ' ' << v[7] << ' ' << v[8] << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const diagTensor& v
    )
    {
        os  << v[0] << ' ' << v[1] << ' ' << v[2] << nl;
    }


    template<>
    inline void Foam::rawSurfaceWriter::writeData
    (
        Ostream& os,
        const tensor& v
    )
    {
        os  << v[0] << ' ' << v[1] << ' ' << v[2] << ' '
            << v[3] << ' ' << v[4] << ' ' << v[5] << ' '
            << v[6] << ' ' << v[7] << ' ' << v[8] << nl;
    }

}


template<class Type>
void Foam::rawSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/fieldName + '_' + surfaceName + ".raw");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    // header
    os  << "# " << fieldName;
    if (isNodeValues)
    {
        os  << "  POINT_DATA ";
    }
    else
    {
        os  << "  FACE_DATA ";
    }

    // header
    writeHeader(os, fieldName, values);

    // values
    if (isNodeValues)
    {
        forAll(values, elemI)
        {
            writeLocation(os, points, elemI);
            writeData(os, values[elemI]);
        }
    }
    else
    {
        forAll(values, elemI)
        {
            writeLocation(os, points, faces, elemI);
            writeData(os, values[elemI]);
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rawSurfaceWriter::rawSurfaceWriter()
:
    surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rawSurfaceWriter::~rawSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rawSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/surfaceName + ".raw");

    if (verbose)
    {
        Info<< "Writing geometry to " << os.name() << endl;
    }


    // header
    os  << "# geometry NO_DATA " << faces.size() << nl
        << "#  x  y  z" << nl;

    // Write faces centres
    forAll(faces, elemI)
    {
        writeLocation(os, points, faces, elemI);
        os  << nl;
    }

    os  << nl;
}


// create write methods
defineSurfaceWriterWriteFields(Foam::rawSurfaceWriter);


// ************************************************************************* //
