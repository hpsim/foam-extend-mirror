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
    faMatrix<Type>

Description
    Finite-Area matrix.

SourceFiles
    faMatrix.C
    faMatrixSolve.C

\*---------------------------------------------------------------------------*/

#ifndef faMatrix_H
#define faMatrix_H

#include "areaFields.H"
#include "edgeFields.H"
#include "lduMatrix.H"
#include "tmp.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class faMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const faMatrix<Type>&);


/*---------------------------------------------------------------------------*\
                           Class faMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faMatrix
:
    public refCount,
    public lduMatrix
{
    // Private data

        // Reference to GeometricField<Type, faPatchField, areaMesh>
        const GeometricField<Type, faPatchField, areaMesh>& psi_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Source term
        Field<Type> source_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //  for internal faces
        FieldField<Field, Type> internalCoeffs_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //  for boundary faces
        FieldField<Field, Type> boundaryCoeffs_;


        //- Face flux field for non-orthogonal correction
        mutable GeometricField<Type, faePatchField, edgeMesh>
            *faceFluxCorrectionPtr_;


    // Private member functions

        //- Add patch contribution to internal field
        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void addToInternalField
        (
            const unallocLabelList& addr,
            const tmp<Field<Type2> >& tpf,
            Field<Type2>& intf
        ) const;

        //- Subtract patch contribution from internal field
        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void subtractFromInternalField
        (
            const unallocLabelList& addr,
            const tmp<Field<Type2> >& tpf,
            Field<Type2>& intf
        ) const;


        // Matrix completion functionality

            void addBoundaryDiag
            (
                scalarField& diag,
                const direction cmpt
            ) const;

            void addCmptAvBoundaryDiag(scalarField& diag) const;

            void addBoundarySource
            (
                Field<Type>& source,
                const bool couples = true
            ) const;


public:

    //- Solver class returned by the solver function
    class faSolver
    {
        faMatrix<Type>& faMat_;

        autoPtr<lduSolver> solver_;

    public:

        // Constructors

            faSolver(faMatrix<Type>& faMat, autoPtr<lduSolver> sol)
            :
                faMat_(faMat),
                solver_(sol)
            {}


        // Member functions

            //- Solve returning the solution statistics.
            //  Solver controls read from dictionary
            lduSolverPerformance solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Solver controls read from faSolution
            lduSolverPerformance solve();
    };


    ClassName("faMatrix");


    // Constructors

        //- Construct given a field to solve for
        faMatrix
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const dimensionSet&
        );

        //- Construct as copy
        faMatrix(const faMatrix<Type>&);

        //- Construct from Istream given field to solve for
        faMatrix
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            Istream&
        );


    //- Destructor
    virtual ~faMatrix();


    // Member functions

        // Access

            const GeometricField<Type, faPatchField, areaMesh>& psi() const
            {
                return psi_;
            }

            const dimensionSet& dimensions() const
            {
                return dimensions_;
            }

            Field<Type>& source()
            {
                return source_;
            }

            const Field<Type>& source() const
            {
                return source_;
            }

            //- faBoundary scalar field containing pseudo-matrix coeffs
            //  for internal cells
            FieldField<Field, Type>& internalCoeffs()
            {
                return internalCoeffs_;
            }

            //- faBoundary scalar field containing pseudo-matrix coeffs
            //  for boundary cells
            FieldField<Field, Type>& boundaryCoeffs()
            {
                return boundaryCoeffs_;
            }


            //- Declare return type of the faceFluxCorrectionPtr() function
            typedef GeometricField<Type, faePatchField, edgeMesh>
                *edgeTypeFieldPtr;

            //- Return pointer to face-flux non-orthogonal correction field
            edgeTypeFieldPtr& faceFluxCorrectionPtr()
            {
                return faceFluxCorrectionPtr_;
            }


        // Operations

            //- Set solution in given cells and eliminate corresponding
            //  equations from the matrix
            void setValues
            (
                const labelList& faces,
                const Field<Type>& values
            );

            //- Set reference level for solution
            void setReference
            (
                const label facei,
                const Type& value
            );

            //- Set reference level for a component of the solution
            //  on a given patch face
            void setComponentReference
            (
                const label patchi,
                const label facei,
                const direction cmpt,
                const scalar value
            );

            //- Relax matrix (for steady-state solution).
            //  alpha = 1 : diagonally equal
            //  alpha < 1 :    ,,      dominant
            //  alpha = 0 : do nothing
            void relax(const scalar alpha);

            //- Relax matrix (for steadty-state solution).
            //  alpha is read from controlDict
            void relax();

            //- Solve returning the solution statistics.
            //  Solver controls read Istream
            lduSolverPerformance solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Solver controls read from faSolution
            lduSolverPerformance solve();

            //- Return the matrix residual
            tmp<Field<Type> > residual() const;

            //- Return the matrix diagonal
            tmp<scalarField> D() const;

            //- Return the central coefficient
            tmp<areaScalarField> A() const;

            //- Return the H operation source
            tmp<GeometricField<Type, faPatchField, areaMesh> > H() const;

            //- Return the face-flux field from the matrix
            tmp<GeometricField<Type, faePatchField, edgeMesh> > flux() const;


    // Member operators

        void operator=(const faMatrix<Type>&);
        void operator=(const tmp<faMatrix<Type> >&);

        void negate();

        void operator+=(const faMatrix<Type>&);
        void operator+=(const tmp<faMatrix<Type> >&);

        void operator-=(const faMatrix<Type>&);
        void operator-=(const tmp<faMatrix<Type> >&);

        void operator+=(const GeometricField<Type,faPatchField,areaMesh>&);
        void operator+=(const tmp<GeometricField<Type,faPatchField,areaMesh> >&);

        void operator-=(const GeometricField<Type,faPatchField,areaMesh>&);
        void operator-=(const tmp<GeometricField<Type,faPatchField,areaMesh> >&);

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const areaScalarField&);
        void operator*=(const tmp<areaScalarField>&);

        void operator*=(const dimensioned<scalar>&);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const faMatrix<Type>&
        );
};


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const faMatrix<Type>&,
    const faMatrix<Type>&,
    const char*
);

template<class Type>
void checkMethod
(
    const faMatrix<Type>&,
    const GeometricField<Type, faPatchField, areaMesh>&,
    const char*
);

template<class Type>
void checkMethod
(
    const faMatrix<Type>&,
    const GeometricField<Type, faPatchField, areaMesh>&,
    const char*
);

template<class Type>
void checkMethod
(
    const faMatrix<Type>&,
    const dimensioned<Type>&,
    const char*
);


//- Solve returning the solution statistics given convergence tolerance
//  Solver controls read Istream
template<class Type>
lduSolverPerformance solve(faMatrix<Type>&, Istream&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Solver controls read Istream
template<class Type>
lduSolverPerformance solve(const tmp<faMatrix<Type> >&, Istream&);


//- Solve returning the solution statistics given convergence tolerance
//  Solver controls read faSolution
template<class Type>
lduSolverPerformance solve(faMatrix<Type>&);


//- Solve returning the solution statistics given convergence tolerance,
//  deleting temporary matrix after solution.
//  Solver controls read faSolution
template<class Type>
lduSolverPerformance solve(const tmp<faMatrix<Type> >&);


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

// Unary operators for negation
template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&
);

// Operators for faMatrix and faMatrix
template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >&,
    const tmp<faMatrix<Type> >&
);


// Operators for faMatrix and GeometricField
template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const GeometricField<Type, faPatchField, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const GeometricField<Type, faPatchField, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const GeometricField<Type, faPatchField, areaMesh>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const GeometricField<Type, faPatchField, areaMesh>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const GeometricField<Type, faPatchField, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const GeometricField<Type, faPatchField, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const GeometricField<Type, faPatchField, areaMesh>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const GeometricField<Type, faPatchField, areaMesh>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
    const tmp<faMatrix<Type> >&
);


// Operators for faMatrix and DimensionedField
template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const DimensionedField<Type, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const DimensionedField<Type, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const tmp<DimensionedField<Type, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const tmp<DimensionedField<Type, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const DimensionedField<Type, areaMesh>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const DimensionedField<Type, areaMesh>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<DimensionedField<Type, areaMesh> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<DimensionedField<Type, areaMesh> >&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const DimensionedField<Type, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const DimensionedField<Type, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const tmp<DimensionedField<Type, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const tmp<DimensionedField<Type, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const DimensionedField<Type, areaMesh>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const DimensionedField<Type, areaMesh>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<DimensionedField<Type, areaMesh> >&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<DimensionedField<Type, areaMesh> >&,
    const tmp<faMatrix<Type> >&
);


// Operators for faMatrix and dimensionedType
template<class Type>
tmp<faMatrix<Type> > operator+
(
    const faMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const faMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const tmp<faMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator+
(
    const dimensioned<Type>&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const tmp<faMatrix<Type> >&,
    const dimensioned<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator-
(
    const dimensioned<Type>&,
    const tmp<faMatrix<Type> >&
);


// Subtraction operators

// Operators for faMatrix and GeometricField
template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>&,
    const GeometricField<Type, faPatchField, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >&,
    const GeometricField<Type, faPatchField, areaMesh>&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>&,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >&,
    const tmp<GeometricField<Type, faPatchField, areaMesh> >&
);

// Operators for faMatrix and dimensionedType
template<class Type>
tmp<faMatrix<Type> > operator==
(
    const faMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator==
(
    const tmp<faMatrix<Type> >&,
    const dimensioned<Type>&
);


// Multiplication operators

// Operators for faMatrix and areaScalarField
template<class Type>
tmp<faMatrix<Type> > operator*
(
    const areaScalarField&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const areaScalarField&,
    const tmp<faMatrix<Type> >&
);

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const tmp<areaScalarField>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const tmp<areaScalarField>&,
    const tmp<faMatrix<Type> >&
);

// Operators for faMatrix and dimensionedScalar
template<class Type>
tmp<faMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const faMatrix<Type>&
);

template<class Type>
tmp<faMatrix<Type> > operator*
(
    const dimensioned<scalar>&,
    const tmp<faMatrix<Type> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
