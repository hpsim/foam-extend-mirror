/*---------------------------------------------------------------------------*\
  ========Merkle=                 |
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
    Foam::phaseChangeTwoPhaseMixtures::Merkle

Description
    Merkle cavitation model.

    Reference:
    @verbatim
        C. L. Merkle, J. Feng, and P. E. O. Buelow,
        "Computational modeling of the dynamics of sheet cavitation",
        in Proceedings Third International Symposium on Cavitation
        Grenoble, France 1998.
   @verbatim

SourceFiles
    Merkle.C

\*--------------------------------------------------------------------*/

#ifndef Merkle_H
#define Merkle_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Merkle
\*--------------------------------------------------------------------*/

class Merkle
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
    TypeName("Merkle");


    // Constructors

        //- construct from components
        Merkle
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    // Destructor

        virtual ~Merkle()
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

        //- Correct the Merkle phaseChange model
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
