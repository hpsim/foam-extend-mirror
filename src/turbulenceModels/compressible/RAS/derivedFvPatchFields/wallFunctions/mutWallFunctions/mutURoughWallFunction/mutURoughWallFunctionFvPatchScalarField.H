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
    Foam::compressible::RASModels::mutURoughWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent (dynamic) viscosity condition
    when using wall functions for rough walls, based on velocity.

Usage
    \table
        Property     | Description             | Required    | Default value
        roughnessHeight | roughness height     | yes         |
        roughnessConstant | roughness constanr | yes         |
        roughnessFactor | scaling factor       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mutURoughWallFunction;
        roughnessHeight 1e-5;
        roughnessConstant 0.5;
        roughnessFactor 1;
    }
    \endverbatim

See also
    Foam::mutWallFunctionFvPatchScalarField

SourceFiles
    mutURoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mutURoughWallFunctionFvPatchScalarField_H
#define mutURoughWallFunctionFvPatchScalarField_H

#include "mutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class mutURoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutURoughWallFunctionFvPatchScalarField
:
    public mutWallFunctionFvPatchScalarField
{
    // Private Data

        // Roughness model parameters

            //- Height
            scalar roughnessHeight_;

            //- Constant
            scalar roughnessConstant_;

            //- Scale factor
            scalar roughnessFactor_;


    // Protected Member Functions

        //- Calculate yPLus
        virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcMut() const;


public:

    //- Runtime type information
    TypeName("mutURoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutURoughWallFunctionFvPatchScalarField
        //  onto a new patch
        mutURoughWallFunctionFvPatchScalarField
        (
            const mutURoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutURoughWallFunctionFvPatchScalarField
        (
            const mutURoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutURoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mutURoughWallFunctionFvPatchScalarField
        (
            const mutURoughWallFunctionFvPatchScalarField&,
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
                new mutURoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the roughness height scale
            const scalar& roughnessHeight() const
            {
                return roughnessHeight_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalar& roughnessHeight()
            {
                return roughnessHeight_;
            }

            //- Return the roughness constant scale
            const scalar& roughnessConstant() const
            {
                return roughnessConstant_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalar& roughnessConstant()
            {
                return roughnessConstant_;
            }

            //- Return the roughness scale factor
            const scalar& roughnessFudgeFactor() const
            {
                return roughnessFactor_;
            }

            //- Return reference to the roughness scale factor to allow
            //  adjustment
            scalar& roughnessFactor()
            {
                return roughnessFactor_;
            }


        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
