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
    Foam::NonEquilibriumReversibleReaction

Description
    Simple extension of Reaction to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    NonEquilibriumReversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef NonEquilibriumReversibleReaction_H
#define NonEquilibriumReversibleReaction_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class NonEquilibriumReversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ReactionRate>
class NonEquilibriumReversibleReaction
:
    public Reaction<ReactionThermo>
{
    // Private data

        ReactionRate fk_;
        ReactionRate rk_;


public:

    //- Runtime type information
    TypeName("nonEquilibriumReversible");


    // Constructors

        //- Construct from components
        NonEquilibriumReversibleReaction
        (
            const Reaction<ReactionThermo>& reaction,
            const ReactionRate& forwardReactionRate,
            const ReactionRate& reverseReactionRate
        );

        //- Construct as copy given new speciesTable
        NonEquilibriumReversibleReaction
        (
            const NonEquilibriumReversibleReaction
                <ReactionThermo, ReactionRate>&,
            const speciesTable& species
        );

        //- Construct from Istream
        NonEquilibriumReversibleReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            Istream& is
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ReactionThermo> > clone() const
        {
            return autoPtr<Reaction<ReactionThermo> >
            (
                new NonEquilibriumReversibleReaction
                    <ReactionThermo, ReactionRate>(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ReactionThermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction<ReactionThermo> >
            (
                new NonEquilibriumReversibleReaction
                    <ReactionThermo, ReactionRate>
                    (*this, species)
            );
        }


    // Destructor

        virtual ~NonEquilibriumReversibleReaction()
        {}


    // Member Functions

        // NonEquilibriumReversibleReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given formard rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NonEquilibriumReversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
