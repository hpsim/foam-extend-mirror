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

#include "linearSpringDamper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearSpringDamper, 0);
    addToRunTimeSelectionTable
    (
        translationalRestraint,
        linearSpringDamper,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearSpringDamper::linearSpringDamper
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    translationalRestraint(name, dict, sixDOF),
    linSpringCoeffs_(dict.lookup("linearSpring")),
    linDampingCoeffs_(dict.lookup("linearDamping"))
{}


Foam::autoPtr<Foam::translationalRestraint>
Foam::linearSpringDamper::clone() const
{
    return autoPtr<translationalRestraint>
    (
        new linearSpringDamper(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearSpringDamper::~linearSpringDamper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::linearSpringDamper::restrainingForce
(
    const scalar,
    const tensor&,
    const vector& x,
    const vector& u
) const
{
    return - (linSpringCoeffs_ & x) - (linDampingCoeffs_ & u);
}


void Foam::linearSpringDamper::write(Ostream& os) const
{
    os.writeKeyword("type") << tab << type()
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("linearSpring") << tab << linSpringCoeffs_
       << token::END_STATEMENT << nl;
    os.writeKeyword("linearDamping") << tab << linDampingCoeffs_
       << token::END_STATEMENT << endl;
}


// ************************************************************************* //
