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
    Foam::XiGModels::KTS

Description
    Simple Kolmogorov time-scale (KTS) model for the flame-wrinling generation
    rate.

SourceFiles
    KTS.C

\*---------------------------------------------------------------------------*/

#ifndef KTS_H
#define KTS_H

#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                          Class KTS Declaration
\*---------------------------------------------------------------------------*/

class KTS
:
    public XiGModel
{
    // Private data

        scalar GEtaCoef;


    // Private Member Functions

        //- Disallow copy construct
        KTS(const KTS&);

        //- Disallow default bitwise assignment
        void operator=(const KTS&);


public:

    //- Runtime type information
    TypeName("KTS");


    // Constructors

        //- Construct from components
        KTS
        (
            const dictionary& XiGProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );


    // Destructor

        virtual ~KTS();


    // Member Functions

        //- Return the flame-wrinking generation rate
        virtual tmp<volScalarField> G() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiGProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiGModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
