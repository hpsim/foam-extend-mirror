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

#include "translationalRestraint.H"
#include "sixDOFODE.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(translationalRestraint, 0);
    defineRunTimeSelectionTable(translationalRestraint, word);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translationalRestraint::translationalRestraint
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    name_(name),
    sixDOF_(sixDOF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::translationalRestraint::~translationalRestraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::translationalRestraint::name() const
{
    return name_;
}


const Foam::sixDOFODE& Foam::translationalRestraint::sixDOF() const
{
    return sixDOF_;
}


Foam::autoPtr<Foam::translationalRestraint> Foam::translationalRestraint::New
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
{
    const word restraintType(dict.lookup("type"));

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(restraintType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "translationalRestraint::New"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict"
            "\n)"
        )   << "Unknown translation restraint type: " << restraintType
            << endl << endl
            << "Valid translation restraint types are: " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<translationalRestraint>
    (
        cstrIter()
        (
            name,
            dict,
            sixDOF
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const translationalRestraint& tr
)
{
    os << tr.name_ << nl << token::BEGIN_BLOCK << nl;

    tr.write(os);

    os << token::END_BLOCK << endl;

    os.check("Ostream& operator<<(Ostream&, const translationalRestraint&");

    return os;
}


// ************************************************************************* //
