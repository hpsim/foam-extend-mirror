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

\*---------------------------------------------------------------------------*/

#include "fixedValuePointPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatch${FieldType},
    ${typeName}FixedValuePointPatch${FieldType}
);


const char* const ${typeName}FixedValuePointPatch${FieldType}::SHA1sum =
    "${SHA1sum}";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FixedValuePointPatch${FieldType}::
${typeName}FixedValuePointPatch${FieldType}
(
    const pointPatch& p,
    const DimensionedField<${TemplateType}, pointMesh>& iF
)
:
    fixedValuePointPatchField<${TemplateType}>(p, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField\n";
    }
}


${typeName}FixedValuePointPatch${FieldType}::
${typeName}FixedValuePointPatch${FieldType}
(
    const ${typeName}FixedValuePointPatch${FieldType}& ptf,
    const pointPatch& p,
    const DimensionedField<${TemplateType}, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<${TemplateType}>(ptf, p, iF, mapper)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField/mapper\n";
    }
}


${typeName}FixedValuePointPatch${FieldType}::
${typeName}FixedValuePointPatch${FieldType}
(
    const pointPatch& p,
    const DimensionedField<${TemplateType}, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValuePointPatchField<${TemplateType}>(p, iF, dict, valueRequired)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/dictionary\n";
    }
}


${typeName}FixedValuePointPatch${FieldType}::
${typeName}FixedValuePointPatch${FieldType}
(
    const ${typeName}FixedValuePointPatch${FieldType}& ptf
)
:
    fixedValuePointPatchField<${TemplateType}>(ptf)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " as copy\n";
    }
}


${typeName}FixedValuePointPatch${FieldType}::
${typeName}FixedValuePointPatch${FieldType}
(
    const ${typeName}FixedValuePointPatch${FieldType}& ptf,
    const DimensionedField<${TemplateType}, pointMesh>& iF
)
:
    fixedValuePointPatchField<${TemplateType}>(ptf, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum} "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FixedValuePointPatch${FieldType}::
~${typeName}FixedValuePointPatch${FieldType}()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}FixedValuePointPatch${FieldType}::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (${verbose:-false})
    {
        Info<<"updateCoeffs ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${code}
//}}} end code

    this->fixedValuePointPatchField<${TemplateType}>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
