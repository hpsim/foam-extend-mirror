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

Namespace
    fac

Description
    Calculate the divergence of the given field.

SourceFiles
    facDiv.C

\*---------------------------------------------------------------------------*/


#ifndef facNDiv_H
#define facNDiv_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const GeometricField<Type, faePatchField, edgeMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const tmp<GeometricField<Type, faePatchField, edgeMesh> >&
    );


    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word& name
    );


    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const edgeScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const tmp<edgeScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word& name
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const edgeScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > ndiv
    (
        const tmp<edgeScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "facNDiv.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
