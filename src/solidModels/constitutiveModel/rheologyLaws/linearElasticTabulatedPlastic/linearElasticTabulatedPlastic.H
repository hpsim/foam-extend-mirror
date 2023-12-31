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
    linearElasticTabulatedPlastic

Description
    Elasto-plastic material law, with linear elasticity and post-yield
    stress strain is given as tabulated data.

SourceFiles
    linearElasticTabulatedPlastic.C

Author
    Philip Cardiff UCD
    philip.cardiff@gmail.com

\*---------------------------------------------------------------------------*/

#ifndef linearElasticTabulatedPlastic_H
#define linearElasticTabulatedPlastic_H

#include "rheologyLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearElasticTabulatedPlastic
:
    public rheologyLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Modulus of elasticity
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Table of post-yield stress versus plastic strain
        // For large strain total Lagrangian models, it expects
        // the 2nd Piola-Kirchhoff stress and Green strain.
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- delta for numerical differentiation used to calculate
        //- local plastic modulus by the finite difference method
        scalar numDiffDelta_;

        //- Initial Yield stress
        dimensionedScalar sigmaY_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticTabulatedPlastic(const linearElasticTabulatedPlastic&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticTabulatedPlastic&);


public:

    //- Runtime type information
    TypeName("linearElasticTabulatedPlastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElasticTabulatedPlastic
        (
            const word& name,
            const volSymmTensorField& sigma,
            const dictionary& dict
        );


    // Destructor

        virtual ~linearElasticTabulatedPlastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const;

        //- Return density
        virtual tmp<volScalarField> rho(scalar t) const
        {
            return rho();
        }

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(scalar t) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const
        {
            return nu();
        }

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const
        {
            notImplemented(type() + "::J(scalar t)");

            return 1.0/E(t);
        }

        //- Return yield stress
        virtual tmp<volScalarField> sigmaY() const;

        //- Return yield stress
        virtual scalar sigmaY
        (const scalar epsilonPEq, const label cellID) const;

        //- If the rheologyLaw needs plasticity procedures
        virtual bool plasticityModelNeeded() const
        {
            return true;
        };

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep() const;

        //- Return plastic modulus as a function of plastic equivalent strain
        virtual tmp<volScalarField> Ep(const volScalarField& epsilonPEq) const;

        //- Correct the rheological model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
