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
    turbulentTemperatureCoupledBaffleFvPatchScalarField

Description
    Harmonic fixed value boundary condition for temperature, to be used
    for heat-transfer on back-to-back baffles.

    If my temperature is T1, heat conductivity Kappa1 and neighbour is T2,Kappa2

    both sides get fixedValue (Kappa1/dx1*T1 + Kappa2/dx2*T2)/(Kappa1/dx1+Kappa2/dx2)

    Example usage:
        myInterfacePatchName
        {
            type                compressible::turbulentTemperatureCoupledBaffle;
            neighbourFieldName  T;
            Kappa                   Kappa; // or none
            value               uniform 300;
        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: if Kappa is "none" looks up RASModel and basicThermo, otherwise expects
    the solver to calculate a 'Kappa' field.

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

    Note: lags interface data so both sides use same data.
    - problem: schedule to calculate average would interfere
    with standard processor swaps.
    - so: updateCoeffs sets both to same Twall. Only need to do
    this for last outer iteration but don't have access to this.

SourceFiles
    turbulentTemperatureCoupledBaffleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCoupledBaffleFvPatchScalarField_H
#define turbulentTemperatureCoupledBaffleFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class turbulentTemperatureCoupledBaffleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCoupledBaffleFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

        //- Name of thermal conductivity field
        const word KappaName_;


    // Private Member Functions

        //- Am I or neighbour owner of interface
        bool interfaceOwner(const polyMesh&, const polyPatch&) const;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureCoupledBaffle");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleFvPatchScalarField onto a new patch
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledBaffleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleFvPatchScalarField&,
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
                new turbulentTemperatureCoupledBaffleFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding Kappa field
        tmp<scalarField> Kappa() const;

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
