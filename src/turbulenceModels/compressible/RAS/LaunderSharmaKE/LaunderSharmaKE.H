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
    Foam::compressible::RASModels::LaunderSharmaKE

Description
    Launder and Sharma low-Reynolds k-epsilon turbulence model for
    compressible and combusting flows.

    The default model coefficients correspond to the following:
    @verbatim
        LaunderSharmaKECoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;
            alphah      1.0;    // only for compressible
            alphahk     1.0;    // only for compressible
            alphaEps    0.76923;
        }
    @endverbatim

SourceFiles
    LaunderSharmaKE.C
    LaunderSharmaKECorrect.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleLaunderSharmaKE_H
#define compressibleLaunderSharmaKE_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class LaunderSharmaKE Declaration
\*---------------------------------------------------------------------------*/

class LaunderSharmaKE
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar Prt_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField mut_;
            volScalarField alphat_;


    // Private member functions

        tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;


public:

    //- Runtime type information
    TypeName("LaunderSharmaKE");

    // Constructors

        //- Construct from components
        LaunderSharmaKE
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~LaunderSharmaKE()
    {}


    // Member Functions

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", mut_/sigmak_ + mu())
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
        virtual tmp<volSymmTensorField> R() const;

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
