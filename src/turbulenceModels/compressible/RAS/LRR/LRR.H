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
    Foam::compressible::RASModels::LRR

Description
    Launder, Reece and Rodi Reynolds-stress turbulence model for
    compressible flows.

    The default model coefficients correspond to the following:
    @verbatim
        LRRCoeffs
        {
            Cmu         0.09;
            Clrr1       1.8;
            Clrr2       0.6;
            C1          1.44;
            C2          1.92;
            Cs          0.25;
            Ceps        0.15;
            Prt         1.0;        // only for compressible
            sigmaEps    1.3;
            sigmaR      0.81967;    // only for compressible
            couplingFactor  0.0;    // only for incompressible
        }
    @endverbatim

SourceFiles
    LRR.C
    LRRcorrect.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleLRR_H
#define compressibleLRR_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LRR Declaration
\*---------------------------------------------------------------------------*/

class LRR
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;

            dimensionedScalar Clrr1_;
            dimensionedScalar Clrr2_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;

            dimensionedScalar couplingFactor_;

            dimensionedScalar sigmaR_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar Prt_;


        // Fields

            volSymmTensorField R_;
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField mut_;
            volScalarField alphat_;


public:

    //- Runtime type information
    TypeName("LRR");

    // Constructors

        //- Construct from components
        LRR
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~LRR()
    {}


    // Member Functions

        //- Return the effective diffusivity for R
        tmp<volScalarField> DREff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DREff", mut_/sigmaR_ + mu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", mut_/sigmaEps_ + mu())
            );
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the turbulence thermal diffusivity
        virtual tmp<volScalarField> alphat() const
        {
            return alphat_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const
        {
            return R_;
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
