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

Description
    Template for use with dynamic code generation of a
    fixedValue pointPatchField.

    - without state

SourceFiles
    fixedValuePointPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchTemplate${FieldType}_H
#define fixedValuePointPatchTemplate${FieldType}_H

#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                A templated FixedValuePointPatchField
\*---------------------------------------------------------------------------*/

class ${typeName}FixedValuePointPatch${FieldType}
:
    public fixedValuePointPatchField<${TemplateType}>
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&,
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&
        );

        //- Construct and return a clone
        virtual autoPtr< pointPatchField<${TemplateType}> > clone() const
        {
            return autoPtr< pointPatchField<${TemplateType}> >
            (
                new ${typeName}FixedValuePointPatch${FieldType}(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr< pointPatchField<${TemplateType}> > clone
        (
            const DimensionedField<${TemplateType}, pointMesh>& iF
        ) const
        {
            return autoPtr< pointPatchField<${TemplateType}> >
            (
                new ${typeName}FixedValuePointPatch${FieldType}(*this, iF)
            );
        }


    //- Destructor
    virtual ~${typeName}FixedValuePointPatch${FieldType}();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
