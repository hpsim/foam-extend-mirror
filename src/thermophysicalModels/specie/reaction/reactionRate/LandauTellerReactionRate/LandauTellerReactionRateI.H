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
inline LandauTellerReactionRate::LandauTellerReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const scalar B,
    const scalar C
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    B_(B),
    C_(C)
{}


//- Construct from Istream
inline LandauTellerReactionRate::LandauTellerReactionRate
(
    const speciesTable&,
    Istream& is
)
:
    A_(readScalar(is.readBegin("LandauTellerReactionRate(Istream&)"))),
    beta_(readScalar(is)),
    Ta_(readScalar(is)),
    B_(readScalar(is)),
    C_(readScalar(is))
{
    is.readEnd("LandauTellerReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar LandauTellerReactionRate::operator()
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

    if (mag(B_) > VSMALL)
    {
        expArg += B_/cbrt(T);
    }

    if (mag(C_) > VSMALL)
    {
        expArg += C_/pow(T, 2.0/3.0);
    }

    if (mag(expArg) > VSMALL)
    {
        lta *= exp(expArg);
    }

    return lta;
}


inline Ostream& operator<<(Ostream& os, const LandauTellerReactionRate& arr)
{
    os  << token::BEGIN_LIST
        << arr.A_ << token::SPACE << arr.beta_ << token::SPACE << arr.Ta_
        << token::SPACE << arr.B_ << token::SPACE << arr.C_
        << token::END_LIST;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
