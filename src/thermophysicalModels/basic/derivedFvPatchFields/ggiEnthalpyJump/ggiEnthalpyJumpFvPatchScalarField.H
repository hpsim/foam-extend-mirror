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
    Foam::ggiEnthalpyJumpFvPatchScalarField

Description
    Foam::ggiEnthalpyJumpFvPatchScalarField

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Gregor Cvijetic, FMENA Zagreb.
    Hrvoje Jasak, Wikki Ltd.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

SourceFiles
    ggiEnthalpyJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ggiEnthalpyJumpFvPatchScalarField_H
#define ggiEnthalpyJumpFvPatchScalarField_H

#include "jumpGgiFvPatchFields.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ggiEnthalpyJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ggiEnthalpyJumpFvPatchScalarField
:
    public jumpGgiFvPatchScalarField
{
    // Private data

        //- Is the patch on the rotating side?
        Switch rotating_;

        //- "jump" field
        scalarField jump_;


public:

    //- Runtime type information
    TypeName("ggiEnthalpyJump");


// Constructors

        //- Construct from patch and internal field
        ggiEnthalpyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ggiEnthalpyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given ggiEnthalpyJumpFvPatchScalarField
        //  onto a new patch
        ggiEnthalpyJumpFvPatchScalarField
        (
            const ggiEnthalpyJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        ggiEnthalpyJumpFvPatchScalarField
        (
            const ggiEnthalpyJumpFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ggiEnthalpyJumpFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new ggiEnthalpyJumpFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Is the patch rotating
            bool rotating() const
            {
                return rotating_;
            }

            //- Return the "jump" across the patch as a "half" field
            virtual tmp<scalarField> jump() const
            {
                return jump_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
