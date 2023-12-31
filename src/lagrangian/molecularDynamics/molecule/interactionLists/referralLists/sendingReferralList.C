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

\*----------------------------------------------------------------------------*/

#include "sendingReferralList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sendingReferralList::sendingReferralList()
:
    labelList(),
    destinationProc_(-1)
{}


Foam::sendingReferralList::sendingReferralList
(
    const label destinationProc,
    const labelList& cellsToSend
)
:
    labelList(cellsToSend),
    destinationProc_(destinationProc)
{}


Foam::sendingReferralList::sendingReferralList
(
    const sendingReferralList& rL
)
:
    labelList(rL),
    destinationProc_(rL.destinationProc())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sendingReferralList::~sendingReferralList()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::sendingReferralList::operator=(const sendingReferralList& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::distribution::operator=(const Foam::distribution&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    labelList::operator=(rhs);

    destinationProc_ = rhs.destinationProc();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool operator==
(
    const Foam::sendingReferralList& a,
    const Foam::sendingReferralList& b
)
{
    // Trivial reject: lists are different size
    if (a.size() != b.size())
    {
        return false;
    }

    // Or if source processors are not the same.
    if (a.destinationProc() != b.destinationProc())
    {
        return false;
    }

    Foam::List<bool> fnd(a.size(), false);

    forAll (b, bI)
    {
        Foam::label curLabel = b[bI];

        bool found = false;

        forAll (a, aI)
        {
            if (a[aI] == curLabel)
            {
                found = true;
                fnd[aI] = true;
                break;
            }
        }

        if (!found)
        {
            return false;
        }
    }

    // check if all labels on a were marked
    bool result = true;

    forAll (fnd, aI)
    {
        result = (result && fnd[aI]);
    }

    return result;
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    sendingReferralList& sRL
)
{
    is  >> sRL.destinationProc_ >> static_cast<labelList&>(sRL);

    is.check("Istream& operator<<(Istream& f, const sendingReferralList& sRL");

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sendingReferralList& rL
)
{
    os  << rL.destinationProc() << token::SPACE
        << static_cast< const labelList& >(rL);

    os.check("Ostream& operator<<(Ostream& f, const sendingReferralList& rL");

    return os;
}


// ************************************************************************* //
