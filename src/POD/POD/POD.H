/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Namespace
    POD

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Proper Orthogonal Decomposition of fields: helper functions

    Note: projection of volFields can be weighted with V/gSum(V)
    but I cannot prove it makes a difference in practical cases.

SourceFiles
    POD.C
    PODIO.C

\*---------------------------------------------------------------------------*/

#ifndef POD_H
#define POD_H

#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace POD functions Declaration
\*---------------------------------------------------------------------------*/

namespace POD
{
    //- Dot product of snapshots across all fields/components
    template<class Type>
    scalar projection(const Field<Type>& a, const Field<Type>& b);

    template<class Type>
    scalar projection(const tmp<Field<Type> >& ta, const Field<Type>& b);

    template<class Type>
    scalar projection(const Field<Type>& a, const tmp<Field<Type> >& tb);

    template<class Type>
    scalar projection(const tmp<Field<Type> >& ta, const tmp<Field<Type> >& tb);


    //- Dot product of snapshots across all fields/components
    template<class Type>
    scalar projection
    (
        const GeometricField<Type, fvPatchField, volMesh>& a,
        const GeometricField<Type, fvPatchField, volMesh>& b
    );

    template<class Type>
    scalar projection
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& ta,
        const GeometricField<Type, fvPatchField, volMesh>& b
    );

    template<class Type>
    scalar projection
    (
        const GeometricField<Type, fvPatchField, volMesh>& a,
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& tb
    );

    template<class Type>
    scalar projection
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& ta,
        const tmp<GeometricField<Type, fvPatchField, volMesh> >& tb
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "POD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
