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
    Foam::immersedBoundaryFvsPatchField

Description
    Foam::immersedBoundaryFvsPatchField

Author
    Hrvoje Jasak

SourceFiles
    immersedBoundaryFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvsPatchField_H
#define immersedBoundaryFvsPatchField_H

#include "fvsPatchField.H"
#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class immersedBoundaryFvPatchField;

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class immersedBoundaryFvsPatchField
:
    public fvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the processor patch
        const immersedBoundaryFvPatch& ibPatch_;


    // Private Member Functions

        //- Update field size to match the patch
        //  Return true if size has changed
        bool updateSize();


public:

    //- Runtime type information
    TypeName(immersedBoundaryFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        immersedBoundaryFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedBoundaryFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given immersedBoundaryFvsPatchField
        //  onto a new patch
        immersedBoundaryFvsPatchField
        (
            const immersedBoundaryFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        immersedBoundaryFvsPatchField
        (
            const immersedBoundaryFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new immersedBoundaryFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        immersedBoundaryFvsPatchField
        (
            const immersedBoundaryFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new immersedBoundaryFvsPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~immersedBoundaryFvsPatchField()
    {}


    // Member functions

        // Access

            //- Return reference to immersed boundary patch
            const immersedBoundaryFvPatch& ibPatch() const
            {
                return ibPatch_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvsPatchField onto this fvsPatchField
            virtual void rmap
            (
                const fvsPatchField<Type>&,
                const labelList&
            );

            //- Update on mesh motion
            virtual void updateOnMotion();


        // Evaluation functions

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::blocking
            );


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        // Note: All assignment operators can check size and resize
        // Other operators cannot, as they assume an existing valid value
        // HJ, 13/Dec/2017


        virtual void operator=(const UList<Type>&);
        virtual void operator=(const fvsPatchField<Type>&);
        virtual void operator=(const fvPatchField<Type>&);

        virtual void operator=(const Type&);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const fvsPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
