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
    Foam::cyclicFaPatchField

Description
    Foam::cyclicFaPatchField

SourceFiles
    cyclicFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFaPatchField_H
#define cyclicFaPatchField_H

#include "coupledFaPatchField.H"
#include "cyclicLduInterfaceField.H"
#include "cyclicFaPatch.H"
#include "areaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicFaPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicFaPatchField
:
    virtual public cyclicLduInterfaceField,
    public coupledFaPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicFaPatch& cyclicPatch_;


    // Private member functions

        //- Return neighbour side field given internal fields
        template<class Type2>
        tmp<Field<Type2> > neighbourSideField
        (
            const Field<Type2>&
        ) const;


public:

    //- Runtime type information
    TypeName(cyclicFaPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicFaPatchField onto a new patch
        cyclicFaPatchField
        (
            const cyclicFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        cyclicFaPatchField
        (
            const cyclicFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >
            (
                new cyclicFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicFaPatchField
        (
            const cyclicFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type> >
            (
                new cyclicFaPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Retirn local reference cast into the cyclic patch
            const cyclicFaPatch& cyclicPatch() const
            {
                return cyclicPatch_;
            }


        // Evaluation functions

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type> > patchNeighbourField() const;


        // Coupled interface functionality

            //- Transform neighbour field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const
            {
                // See comments in cyclicLduInterfaceField
                // HJ, 3/Sep/2015
//                 cyclicLduInterfaceField::transformCoupleField(f, cmpt);
                if (doTransform())
                {
                    label sizeby2 = f.size()/2;

                    Type powFwdTransform = transform
                    (
                        this->forwardT()[0],
                        pTraits<Type>::one
                    );

                    scalar forwardScale =
                        pow(component(powFwdTransform, cmpt), rank());

                    Type powRevTransform = transform
                    (
                        this->reverseT()[0],
                        pTraits<Type>::one
                    );

                    scalar reverseScale =
                        pow(component(powRevTransform, cmpt), rank());

                    for (label facei = 0; facei < sizeby2; facei++)
                    {
                        f[facei] *= forwardScale;
                        f[facei + sizeby2] *= reverseScale;
                    }
                }
            }

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        // Block coupled interface functionality

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field<Type>&,
                Field<Type>&,
                const BlockLduMatrix<Type>&,
                const CoeffField<Type>&,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const
            {
                notImplemented
                (
                    "cyclicFaPatchField<Type>::"
                    "updateInterfaceMatrix for block matrices"
                );
            }


        //- Cyclic coupled interface functions

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return !(cyclicPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cyclicFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
