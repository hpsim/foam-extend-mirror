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
    Foam::phaseChangeTwoPhaseMixtures::Kunz

Description
    Kunz cavitation model slightly modified so that the condensation term
    is switched off when the pressure is less than the saturation vapour
    pressure.  This change allows the condensation term to be formulated as
    a coefficient multiplying (p - p_sat) so that it can be included as an
    implicit term in the pressure equation.

    Reference:
    @verbatim
        Kunz, R.F., Boger, D.A., Stinebring, D.R., Chyczewski, Lindau. J.W.,
        Gibeling, H.J., Venkateswaran, S., Govindan, T.R.,
        “A Preconditioned Implicit Method for Two-Phase Flows with Application
         to Cavitation Prediction,”
        Computers and Fluids,
        29(8):849-875, 2000.
    @verbatim

SourceFiles
    Kunz.C

\*--------------------------------------------------------------------*/

#ifndef Kunz_H
#define Kunz_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Kunz
\*--------------------------------------------------------------------*/

class Kunz
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

        dimensionedScalar UInf_;
        dimensionedScalar tInf_;
        dimensionedScalar Cc_;
        dimensionedScalar Cv_;

        dimensionedScalar p0_;

        dimensionedScalar mcCoeff_;
        dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("Kunz");


    // Constructors

        //- construct from components
        Kunz
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        virtual ~Kunz()
        {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField> > mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as an
        //  explicit term for the condensation rate and a coefficient to
        //  multiply (p - pSat) for the vaporisation rate
        virtual Pair<tmp<volScalarField> > mDotP() const;

        //- Correct the Kunz phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
