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

Description
    Static functions in Proper Orthogonal Decomposition

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace POD
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Note: Consider weighting projection with cell volumes

template<class Type>
scalar projection(const Field<Type>& a, const Field<Type>& b)
{
    return cmptSum(gSumCmptProd(a, b));
}


template<class Type>
scalar projection(const tmp<Field<Type> >& ta, const Field<Type>& b)
{
    scalar p = projection(ta(), b);
    ta.clear();

    return p;
}


template<class Type>
scalar projection(const Field<Type>& a, const tmp<Field<Type> >& tb)
{
    scalar p = projection(a, tb());
    tb.clear();

    return p;
}


template<class Type>
scalar projection(const tmp<Field<Type> >& ta, const tmp<Field<Type> >& tb)
{
    scalar p = projection(ta(), tb());
    ta.clear();
    tb.clear();

    return p;
}


template<class Type>
scalar projection
(
    const GeometricField<Type, fvPatchField, volMesh>& a,
    const GeometricField<Type, fvPatchField, volMesh>& b
)
{
    scalar p = projection(a.internalField(), b.internalField());

    // Experimental; reconsider.  HJ, 10/Feb/2021
    // forAll (a.boundaryField(), patchI)
    // {
    //     p += projection(a.boundaryField()[patchI], b.boundaryField()[patchI]);
    // }

    return p;
}


template<class Type>
scalar projection
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& ta,
    const GeometricField<Type, fvPatchField, volMesh>& b
)
{
    const GeometricField<Type, fvPatchField, volMesh>& a = ta();

    scalar p = projection(a, b);
    ta.clear();

    return p;
}


template<class Type>
scalar projection
(
    const GeometricField<Type, fvPatchField, volMesh>& a,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tb
)
{
    const GeometricField<Type, fvPatchField, volMesh>& b = tb();

    scalar p = projection(a, b);
    tb.clear();

    return p;
}


template<class Type>
scalar projection
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& ta,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tb
)
{
    const GeometricField<Type, fvPatchField, volMesh>& a = ta();
    const GeometricField<Type, fvPatchField, volMesh>& b = tb();

    scalar p = projection(a, b);
    ta.clear();
    tb.clear();

    return p;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace POD

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
