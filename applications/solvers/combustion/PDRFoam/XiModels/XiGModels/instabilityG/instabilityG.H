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
    Foam::XiGModels::instabilityG

Description
    Flame-surface instabilityG flame-wrinking generation rate coefficient model
    used in \link XiModel.H \endlink.

    See Technical Report SH/RE/01R for details on the PDR modelling.

SourceFiles
    instabilityG.C

\*---------------------------------------------------------------------------*/

#ifndef instabilityG_H
#define instabilityG_H

#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                          Class instabilityG Declaration
\*---------------------------------------------------------------------------*/

class instabilityG
:
    public XiGModel
{
    // Private data

        //- Flame instabilityG wrinling generation rate coefficient
        dimensionedScalar GIn;

        //- InstabilityG length-scale
        dimensionedScalar lambdaIn;

        //- Xi generation rate model due to all other processes
        autoPtr<XiGModel> XiGModel_;


    // Private Member Functions

        //- Disallow copy construct
        instabilityG(const instabilityG&);

        //- Disallow default bitwise assignment
        void operator=(const instabilityG&);


public:

    //- Runtime type information
    TypeName("instabilityG");


    // Constructors

        //- Construct from components
        instabilityG
        (
            const dictionary& XiGProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~instabilityG();


    // Member Functions

        //- Return the flame-wrinking generation rate
        virtual tmp<volScalarField> G() const;

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiGModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
