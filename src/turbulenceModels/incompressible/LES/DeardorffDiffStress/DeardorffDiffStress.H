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
    Foam::incompressible::LESModels::DeardorffDiffStress

Description
    Differential SGS Stress Equation Model for incompressible flows

    The DSEM uses a model version of the full balance equation for the SGS
    stress tensor to simulate the behaviour of B.
    Thus,
    @verbatim
        d/dt(B) + div(U*B) - div(nuSgs*grad(B))
        =
        P - c1*epsilon/k*B - 0.667*(1 - c1)*epsilon*I - c2*(P - 0.333*trP*I)

    where

        k = 0.5*tr(B),
        epsilon = ce*k^3/2/delta,
        epsilon/k = ce*k^1/2/delta
        P = -(B'L + L'B)
        nuSgs = ck*sqrt(k)*delta
        nuEff = nuSgs + nu
    @endverbatim

SourceFiles
    DeardorffDiffStress.C

\*---------------------------------------------------------------------------*/

#ifndef DeardorffDiffStress_H
#define DeardorffDiffStress_H

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class DeardorffDiffStress Declaration
\*---------------------------------------------------------------------------*/

class DeardorffDiffStress
:
    public GenSGSStress
{
    // Private data

        dimensionedScalar ck_;
        dimensionedScalar cm_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volScalarField& K);

        // Disallow default bitwise copy construct and assignment
        DeardorffDiffStress(const DeardorffDiffStress&);
        DeardorffDiffStress& operator=(const DeardorffDiffStress&);


public:

    //- Runtime type information
    TypeName("DeardorffDiffStress");

    // Constructors

        //- Construct from components
        DeardorffDiffStress
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~DeardorffDiffStress()
    {}


    // Member Functions

        //- Return the effective diffusivity for B
        tmp<volScalarField> DBEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DBEff", nuSgs_ + nu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
