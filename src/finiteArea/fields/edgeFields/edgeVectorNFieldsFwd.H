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
    edgeVectorNFields

Description

SourceFiles
    edgeVectorNFields.C

\*---------------------------------------------------------------------------*/

#ifndef edgeVectorNFieldsFwd_H
#define edgeVectorNFieldsFwd_H

#include "edgeFieldsFwd.H"

#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class edgeMesh;

template<class Type>
class faePatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

#define doMakeTypedefs(type, Type, args...)         \
    typedef GeometricField<type, faePatchField, edgeMesh> edge##Type##Field;

forAllVectorNTypes(doMakeTypedefs)

forAllTensorNTypes(doMakeTypedefs)

forAllDiagTensorNTypes(doMakeTypedefs)

forAllSphericalTensorNTypes(doMakeTypedefs)

#undef doMakeTypedefs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
