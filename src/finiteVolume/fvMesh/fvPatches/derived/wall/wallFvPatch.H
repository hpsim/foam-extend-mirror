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
    Foam::wallFvPatch

Description
    Foam::wallFvPatch

SourceFiles
    wallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wallFvPatch_H
#define wallFvPatch_H

#include "fvPatch.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class wallFvPatch
:
    public fvPatch
{
    // Private Data

        //- Reference to polyPatch
        const wallPolyPatch& wallPolyPatch_;


protected:

    // Protected Member Functions

        // Geometry correction functions

            //- Make patch face non-orthogonality correction vectors
            virtual void makeCorrVecs(fvsPatchVectorField&) const;

            //- Update mesh motion fluxes
            virtual void updatePhi
            (
                DimensionedField<scalar, volMesh>& V,
                DimensionedField<scalar, volMesh>& V0,
                surfaceScalarField& phi
            ) const;


public:

    //- Runtime type information
    TypeName(wallPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        wallFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            wallPolyPatch_(refCast<const wallPolyPatch>(patch))
        {}


    //- Destructor
    virtual ~wallFvPatch() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
