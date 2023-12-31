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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline JanevReactionRate::JanevReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar b[]
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta)
{
    for (int n=0; n<nb_; n++)
    {
        b_[n] = b[n];
    }
}


//- Construct from Istream
inline JanevReactionRate::JanevReactionRate
(
    const speciesTable&,
    Istream& is
)
:
    A_(readScalar(is.readBegin("JanevReactionRate(Istream&)"))),
    beta_(readScalar(is)),
    Ta_(readScalar(is))
{
    for (int n=0; n<nb_; n++)
    {
        is >> b_[n];
    }

    is.readEnd("JanevReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar JanevReactionRate::operator()
(
    const scalar T,
    const scalar,
    const scalarField&
) const
{
    scalar lta = A_;

    if (mag(beta_) > VSMALL)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0.0;

    if (mag(Ta_) > VSMALL)
    {
        expArg -= Ta_/T;
    }

    scalar lnT = log(T);

    for (int n=0; n<nb_; n++)
    {
        expArg += b_[n]*pow(lnT, n);
    }

    lta *= exp(expArg);

    return lta;
}


inline Ostream& operator<<(Ostream& os, const JanevReactionRate& jrr)
{
    os  << token::BEGIN_LIST
        << jrr.A_ << token::SPACE << jrr.beta_ << token::SPACE << jrr.Ta_;

    for (int n=0; n<JanevReactionRate::nb_; n++)
    {
        os  << token::SPACE << jrr.b_[n];
    }

    os << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
