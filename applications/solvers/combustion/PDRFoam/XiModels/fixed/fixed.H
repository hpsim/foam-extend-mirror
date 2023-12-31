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
    Foam::XiModels::fixed

Description
    Fixed value model for Xi. See \link XiModel.H \endlink for more details
    on flame wrinkling modelling.

SourceFiles
    fixed.C

\*---------------------------------------------------------------------------*/

#ifndef fixed_H
#define fixed_H

#include "XiModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class fixed Declaration
\*---------------------------------------------------------------------------*/

class fixed
:
    public XiModel
{
    // Private Member Functions

        //- Disallow copy construct
        fixed(const fixed&);

        //- Disallow default bitwise assignment
        void operator=(const fixed&);


public:

    //- Runtime type information
    TypeName("fixed");


    // Constructors

        //- Construct from components
        fixed
        (
            const dictionary& XiProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~fixed();


    // Member Functions

        //- Correct the flame-wrinking Xi
        virtual void correct()
        {}

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
