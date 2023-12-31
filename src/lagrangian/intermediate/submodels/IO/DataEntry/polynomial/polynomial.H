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
    Foam::polynomial

Description
    Polynomial container data entry for scalars. Items are stored in a list of
    Tuple2's. Data is input in the form, e.g. for an entry <entryName> that
    describes y = x^2 + 2x^3

    @verbatim
        <entryName>   polynomial
        (
            (1    2)
            (2    3)
        );
    @endverbatim

SourceFiles
    polynomial.C

\*---------------------------------------------------------------------------*/

#ifndef polynomial_H
#define polynomial_H

#include "DataEntry.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polynomial;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const polynomial&
);

/*---------------------------------------------------------------------------*\
                        Class polynomial Declaration
\*---------------------------------------------------------------------------*/

class polynomial
:
    public DataEntry<scalar>
{
    // Private data

        //- Polynomial coefficients - list of prefactor, exponent
        List<Tuple2<scalar, scalar> > coeffs_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const polynomial&);


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from entry name and Istream
        polynomial(const word& entryName, Istream& is);

        //- Copy constructor
        polynomial(const polynomial& poly);


    //- Destructor
    virtual ~polynomial();


    // Member Functions

        //- Return polynomial value
        scalar value(const scalar x) const;

        //- Integrate between two (scalar) values
        scalar integrate(const scalar x1, const scalar x2) const;


    //- Ostream Operator
    friend Ostream& operator<<
    (
        Ostream&,
        const polynomial&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
