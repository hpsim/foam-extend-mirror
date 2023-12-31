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
    Foam::referredMolecule

Description

SourceFiles
    referredMoleculeI.H
    referredMolecule.C

\*---------------------------------------------------------------------------*/

#ifndef referredMolecule_H
#define referredMolecule_H

#include "vector.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class referredMolecule;

inline bool operator==(const referredMolecule& a, const referredMolecule& b);
inline bool operator!=(const referredMolecule& a, const referredMolecule& b);

Istream& operator>>(Istream&, referredMolecule&);
Ostream& operator<<(Ostream&, const referredMolecule&);


/*---------------------------------------------------------------------------*\
                      Class referredMolecule Declaration
\*---------------------------------------------------------------------------*/

class referredMolecule
{
    // Private data

        label id_;

        vector position_;

        List<vector> sitePositions_;


public:

    // Constructors

        //- Construct null
        referredMolecule();

        //- Construct from components
        referredMolecule
        (
            const label id,
            const vector& position,
            const List<vector>& sitePositions
        );


    // Destructor

        virtual ~referredMolecule();


    // Member Functions

        // Access

            inline label id() const;

            inline const vector& position() const;

            inline const List<vector>& sitePositions() const;


    // Friend Operators

        inline friend bool operator==
        (
            const referredMolecule& a,
            const referredMolecule& b
        );

        inline friend bool operator!=
        (
            const referredMolecule& a,
            const referredMolecule& b
        );


    // IOstream Operators

        friend Istream& operator>>
        (
            Istream&,
            referredMolecule&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const referredMolecule&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredMoleculeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
