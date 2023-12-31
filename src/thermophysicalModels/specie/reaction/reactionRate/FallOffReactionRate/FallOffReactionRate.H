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
    Foam::FallOffReactionRate

Description
    General class for handling unimolecular/recombination fall-off reactions.

SourceFiles
    FallOffReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef FallOffReactionRate_H
#define FallOffReactionRate_H

#include "thirdBodyEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ReactionRate, class FallOffFunction>
class FallOffReactionRate;

template<class ReactionRate, class FallOffFunction>
inline Ostream& operator<<
(
    Ostream&,
    const FallOffReactionRate<ReactionRate, FallOffFunction>&
);


/*---------------------------------------------------------------------------*\
              Class FallOffReactionRate Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionRate, class FallOffFunction>
class FallOffReactionRate
{
    // Private data

        ReactionRate k0_;
        ReactionRate kInf_;
        FallOffFunction F_;
        thirdBodyEfficiencies thirdBodyEfficiencies_;


public:

    // Constructors

        //- Construct from components
        inline FallOffReactionRate
        (
            const ReactionRate& k0,
            const ReactionRate& kInf,
            const FallOffFunction& F,
            const thirdBodyEfficiencies& tbes
        );

        //- Construct from Istream
        inline FallOffReactionRate
        (
            const speciesTable& species,
            Istream& is
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return ReactionRate::type() + FallOffFunction::type() + "FallOff";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar p,
            const scalarField& c
        ) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionRate, FallOffFunction>
        (
            Ostream&,
            const FallOffReactionRate<ReactionRate, FallOffFunction>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FallOffReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
