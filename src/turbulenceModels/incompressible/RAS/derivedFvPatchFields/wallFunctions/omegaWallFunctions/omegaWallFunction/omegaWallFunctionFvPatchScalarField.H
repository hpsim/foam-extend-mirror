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
    Foam::incompressible::RASModels::omegaWallFunctionFvPatchScalarField

Description
    Provides a wall function boundary condition/constraint on omega

    Computed value is:

        omega = sqrt(omega_vis^2 + omega_log^2)

    where
        omega_vis = omega in viscous region
        omega_log = omega in logarithmic region

    Model described by Eq.(15) of:
    @verbatim
        Menter, F., Esch, T.
        "Elements of Industrial Heat Transfer Prediction"
        16th Brazilian Congress of Mechanical Engineering (COBEM),
        Nov. 2001
    @endverbatim

SourceFiles
    omegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaWallFunctionFvPatchScalarField_H
#define omegaWallFunctionFvPatchScalarField_H

#include "fixedInternalValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class omegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaWallFunctionFvPatchScalarField
:
    public fixedInternalValueFvPatchField<scalar>
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of turbulence kinetic energy field
        word kName_;

        //- Name of turbulence generation field
        word GName_;

        //- Name of laminar viscosity field
        word nuName_;

        //- Name of turbulent viscosity field
        word nutName_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- beta1 coefficient
        scalar beta1_;


    // Private member functions

        //- Check the type of the patch
        void checkType();


protected:

    // Protected Member Functions

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- Return name of turbulence generation field
        inline const word& GName() const
        {
            return GName_;
        }


public:

    //- Runtime type information
    TypeName("omegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // omegaWallFunctionFvPatchScalarField
        //  onto a new patch
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
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
                new omegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
