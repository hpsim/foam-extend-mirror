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

Class
    Foam::fileFormats::VTKsurfaceFormat

Description
    Provide a means of writing VTK legacy format.
    The output is never sorted by zone.

SourceFiles
    VTKsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef VTKsurfaceFormat_H
#define VTKsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "VTKsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                     Class VTKsurfaceFormat Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class VTKsurfaceFormat
:
    public MeshedSurface<Face>,
    public VTKsurfaceFormatCore
{
    // Private Member Functions

        //- Write header information about number of polygon points
        static void writeHeaderPolygons(Ostream&, const UList<Face>&);

        //- Disallow default bitwise copy construct
        VTKsurfaceFormat(const VTKsurfaceFormat<Face>&);

        //- Disallow default bitwise assignment
        void operator=(const VTKsurfaceFormat<Face>&);


public:

    // Constructors

        //- Construct null
        VTKsurfaceFormat();


    // Destructor

        virtual ~VTKsurfaceFormat()
        {}


    // Member Functions

        // Write

            //- Write surface mesh components by proxy
            static void write
            (
                const fileName&, const MeshedSurfaceProxy<Face>&
            );

            //- Write UnsortedMeshedSurface, the output remains unsorted
            static void write
            (
                const fileName&, const UnsortedMeshedSurface<Face>&
            );


            //- Write object
            virtual void write(Ostream& os) const
            {
                write(os, MeshedSurfaceProxy<Face>(*this));
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "VTKsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
