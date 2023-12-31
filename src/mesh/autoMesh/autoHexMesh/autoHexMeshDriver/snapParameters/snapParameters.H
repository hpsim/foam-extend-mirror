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
    Foam::snapParameters

Description
    Simple container to keep together snap specific information.

SourceFiles
    snapParameters.C

\*---------------------------------------------------------------------------*/

#ifndef snapParameters_H
#define snapParameters_H

#include "dictionary.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class snapParameters Declaration
\*---------------------------------------------------------------------------*/

class snapParameters
{
    // Private data

        const label nSmoothPatch_;

        const scalar snapTol_;

        const label nSmoothDispl_;

        const label nSnap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        snapParameters(const snapParameters&);

        //- Disallow default bitwise assignment
        void operator=(const snapParameters&);


public:

    // Constructors

        //- Construct from dictionary - old syntax
        snapParameters(const dictionary& dict, const label dummy);

        //- Construct from dictionary - new syntax
        snapParameters(const dictionary& dict);


    // Member Functions

        // Access

            //- Number of patch smoothing iterations before finding
            //  correspondence to surface
            label nSmoothPatch() const
            {
                return nSmoothPatch_;
            }

            //- Relative distance for points to be attracted by surface
            //  feature point
            //  or edge. True distance is this factor times local
            //  maximum edge length.
            scalar snapTol() const
            {
                return snapTol_;
            }

            //- Number of mesh displacement smoothing iterations.
            label nSmoothDispl() const
            {
                return nSmoothDispl_;
            }

            //- Maximum number of snapping relaxation iterations. Should stop
            //  before upon reaching a correct mesh.
            label nSnap() const
            {
                return nSnap_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
