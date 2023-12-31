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
    Foam::incompressible::RASModels::
    nutSpalartAllmarasStandardWallFunctionFvPatchScalarField

Description
    This boundary condition is the same as nutUWallFunction. We will simply
    derive from nutUWallFunctionFvPatchScalarField to enable backward
    compatibility.

SourceFiles
    nutSpalartAllmarasStandardWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutSpalartAllmarasStandardWallFunctionFvPatchScalarField_H
#define nutSpalartAllmarasStandardWallFunctionFvPatchScalarField_H

#include "nutUWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
  Class nutSpalartAllmarasStandardWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
:
    public nutUWallFunctionFvPatchScalarField
{
public:

    //- Runtime type information
    TypeName("nutSpalartAllmarasStandardWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
        //  onto a new patch
        nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
        (
            const nutSpalartAllmarasStandardWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
        (
            const nutSpalartAllmarasStandardWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
        (
            const nutSpalartAllmarasStandardWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutSpalartAllmarasStandardWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
