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
    Foam::incompressible::LESModels::scaleSimilarity

Description
    General base class for all scale similarity models
    for incompressible flows.

    Since such models do not work without additional eddy viscosity terms,
    this class must be combined with an eddy viscosity model of some form.

SourceFiles
    scaleSimilarity.C

\*---------------------------------------------------------------------------*/

#ifndef scaleSimilarity_H
#define scaleSimilarity_H

#include "LESModel.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class scaleSimilarity Declaration
\*---------------------------------------------------------------------------*/

class scaleSimilarity
:
    virtual public LESModel
{
    // Private data

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        scaleSimilarity(const scaleSimilarity&);
        scaleSimilarity& operator=(const scaleSimilarity&);


public:

    //- Runtime type information
    TypeName("scaleSimilarity");

    // Constructors

        //- Construct from components
        scaleSimilarity
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~scaleSimilarity()
    {}


    // Member Functions

        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const;

        //- Return the SGS turbulent dissipation.
        virtual tmp<volScalarField> epsilon() const;

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const;

        //- Return the deviatoric part of the divergence of Beff
        //  i.e. the additional term in the filtered NSE.
        virtual tmp<fvVectorMatrix> divDevBeff() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>&);

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
