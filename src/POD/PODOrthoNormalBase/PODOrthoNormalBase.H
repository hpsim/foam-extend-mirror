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

Class
    PODOrthoNormalBase

Description
    Establish POD ortho-normal base and interpolation coefficients give a list
    of fields. Size of ortho-normal base is calculated from the desired
    accuracy, e.g. 0.99-0.99999 (in energy terms)

    Orho-normal base is defined on a list of snapshots.  If additional
    non-defining fields are required, they follow the decomposition and weights
    from defining snapshots

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    PODOrthoNormalBase.C

\*---------------------------------------------------------------------------*/

#ifndef PODOrthoNormalBase_H
#define PODOrthoNormalBase_H

#include "volFields.H"
#include "scalarMatrices.H"
#include "PODEigenBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PODOrthoNormalBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PODOrthoNormalBase
{
public:

    // Public typedefs

        //- Defining snapshots field type
        typedef GeometricField<Type, fvPatchField, volMesh> GeoTypeField;


private:

    // Private data

        //- Defining eigen-base
        const PODEigenBase<Type> eigenBase_;

        //- List of defining ortho-normal fields
        PtrList<GeoTypeField> orthoFields_;

        //- Interpolation coefficients
        scalarRectangularMatrix* interpolationCoeffsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PODOrthoNormalBase(const PODOrthoNormalBase&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const PODOrthoNormalBase&) = delete;


public:

    // Constructors

        //- Construct given list of snapshots and accuracy
        PODOrthoNormalBase
        (
            const PtrList<GeoTypeField>& snapshots,
            const scalar accuracy
        );


    // Destructor
    ~PODOrthoNormalBase();


    // Member Functions

        //- Return eigen-base
        inline const PODEigenBase<Type>& eigenBase() const
        {
            return eigenBase_;
        }

        //- Return base size
        inline label baseSize() const
        {
            return orthoFields_.size();
        }

        //- Return n-th ortho-normal base field
        inline const GeoTypeField& orthoField(const label n) const
        {
            return orthoFields_[n];
        }

        //- Return interpolation coefficients (snapshots x orthoVectors)
        inline const scalarRectangularMatrix& interpolationCoeffs() const
        {
            return *interpolationCoeffsPtr_;
        }

        //- Check ortho-normal base: orthogonality and magnitude of snapshots
        void checkBase() const;

        //- Calculate ortho-normal base and fill into orthoFields
        //  Size of orthoFields defines the required size of decomposition
        //  Function is used when further fields are carried on an existing
        //  ortho-normal base
        template<class GeoField>
        void calcOrthoBase
        (
            const PtrList<GeoField>& snapshots,
            PtrList<GeoField>& orthoFields
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PODOrthoNormalBase.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
