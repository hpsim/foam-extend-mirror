/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "error.H"
#include "blockEdge.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockEdge, 0);
    defineRunTimeSelectionTable(blockEdge, Istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockEdge::blockEdge
(
    const pointField& points,
    const label start,
    const label end
)
:
    points_(points),
    start_(start),
    end_(end)
{}


Foam::blockEdge::blockEdge(const pointField& points, Istream& is)
:
    points_(points),
    start_(readLabel(is)),
    end_(readLabel(is))
{}


Foam::blockEdge::blockEdge(const blockEdge& c)
:
    points_(c.points_),
    start_(c.start_),
    end_(c.end_)
{}


Foam::autoPtr<Foam::blockEdge> Foam::blockEdge::New
(
    const pointField& points,
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction
            << "constructing blockEdge"
            << endl;
    }

    const word edgeType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(edgeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown blockEdge type "
            << edgeType << nl << nl
            << "Valid blockEdge types are" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<blockEdge>(cstrIter()(points, is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::blockEdge::appendEndPoints
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherKnots
)
{
    tmp<pointField> tallKnots(new pointField(otherKnots.size() + 2));
    pointField& allKnots = tallKnots();

    // Add start/end knots
    allKnots[0] = points[start];
    allKnots[otherKnots.size() + 1] = points[end];

    // Add intermediate knots
    forAll(otherKnots, knotI)
    {
        allKnots[knotI + 1] = otherKnots[knotI];
    }

    return tallKnots;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::blockEdge::operator=(const blockEdge&)
{
    NotImplemented;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const blockEdge& p)
{
    os << p.start_ << tab << p.end_ << endl;

    return os;
}


// ************************************************************************* //
