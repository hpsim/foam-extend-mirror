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
    Foam::errorEstimate

Description
    Residual error estimation

SourceFiles
    errorEstimate.C

\*---------------------------------------------------------------------------*/

#ifndef errorEstimate_H
#define errorEstimate_H

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class errorEstimate Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class errorEstimate
:
    public refCount
{
    // Private data

        // Reference to GeometricField<Type, fvPatchField, volMesh>
        const GeometricField<Type, fvPatchField, volMesh>& psi_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Cell residual pointer
        Field<Type> residual_;

        //- Normalisation factor
        scalarField normFactor_;


    // Private Member Functions

        //- Return boundary condition types for the error field
        wordList errorBCTypes() const;

public:

    // Static data members

    ClassName("errorEstimate");


    // Constructors

        //- Construct from components
        errorEstimate
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            const dimensionSet& ds,
            const Field<Type>& res,
            const scalarField& norm
        );

        //- Construct as copy
        errorEstimate(const errorEstimate<Type>&);


    // Destructor

        ~errorEstimate();


    // Member Functions

        // Access

            //- Return field
            const GeometricField<Type, fvPatchField, volMesh>& psi() const
            {
                return psi_;
            }

            //- Return residual dimensions
            const dimensionSet& dimensions() const
            {
                return dimensions_;
            }

        // Raw residual (for calculus)

            Field<Type>& res()
            {
                return residual_;
            }

            const Field<Type>& res() const
            {
                return residual_;
            }


        // Error Estimate

            //- Cell residual (volume intensive)
            tmp<GeometricField<Type, fvPatchField, volMesh> > residual() const;

            //- Normalisation factor
            tmp<volScalarField> normFactor() const;

            //- Error estimate
            tmp<GeometricField<Type, fvPatchField, volMesh> > error() const;


    // Member Operators

        void operator=(const errorEstimate<Type>&);
        void operator=(const tmp<errorEstimate<Type> >&);

        void negate();

        void operator+=(const errorEstimate<Type>&);
        void operator+=(const tmp<errorEstimate<Type> >&);

        void operator-=(const errorEstimate<Type>&);
        void operator-=(const tmp<errorEstimate<Type> >&);

        void operator+=(const GeometricField<Type,fvPatchField,volMesh>&);
        void operator+=(const tmp<GeometricField<Type,fvPatchField,volMesh> >&);

        void operator-=(const GeometricField<Type,fvPatchField,volMesh>&);
        void operator-=(const tmp<GeometricField<Type,fvPatchField,volMesh> >&);

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const volScalarField&);
        void operator*=(const tmp<volScalarField>&);

        void operator*=(const dimensioned<scalar>&);


    // Friend Functions

    // Friend Operators
};


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const errorEstimate<Type>&,
    const errorEstimate<Type>&,
    const char*
);

template<class Type>
void checkMethod
(
    const errorEstimate<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&,
    const char*
);

template<class Type>
void checkMethod
(
    const errorEstimate<Type>&,
    const dimensioned<Type>&,
    const char*
);


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const dimensioned<Type>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const dimensioned<Type>&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >&,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >&,
    const dimensioned<Type>&
);


template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const volScalarField&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const volScalarField&,
    const tmp<errorEstimate<Type> >&
);

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const tmp<volScalarField>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const tmp<volScalarField>&,
    const tmp<errorEstimate<Type> >&
);


template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const dimensioned<scalar>&,
    const errorEstimate<Type>&
);

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const dimensioned<scalar>&,
    const tmp<errorEstimate<Type> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "errorEstimate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
