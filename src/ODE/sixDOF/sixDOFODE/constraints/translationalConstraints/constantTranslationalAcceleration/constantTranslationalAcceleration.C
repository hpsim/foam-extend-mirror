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

#include "constantTranslationalAcceleration.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantTranslationalAcceleration, 0);
    addToRunTimeSelectionTable
    (
        translationalConstraint,
        constantTranslationalAcceleration,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantTranslationalAcceleration::constantTranslationalAcceleration
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    translationalConstraint(name, dict, sixDOF),
    dir_(dict.lookup("constraintDirection")),
    a_(readScalar(dict.lookup("translationalAcceleration")))
{
    // Rescale direction
    if (mag(dir_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::constantTranslationalAcceleration::"
            "constantTranslationalAcceleration"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict"
            "\n)"
        )   << "Zero direction specified. This is not allowed."
            << exit(FatalError);
    }
    else
    {
        dir_ /= mag(dir_);
    }
}


Foam::autoPtr<Foam::translationalConstraint>
Foam::constantTranslationalAcceleration::clone() const
{
    return autoPtr<translationalConstraint>
    (
        new constantTranslationalAcceleration(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantTranslationalAcceleration::~constantTranslationalAcceleration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::constantTranslationalAcceleration::matrixContribution
(
    const scalar,
    const tensor&,
    const vector&,
    const vector&
) const
{
    return dir_;
}


Foam::scalar Foam::constantTranslationalAcceleration::sourceContribution
(
    const scalar,
    const tensor&,
    const vector&,
    const vector&
) const
{
    return a_;
}


void Foam::constantTranslationalAcceleration::write(Ostream& os) const
{
    os.writeKeyword("type") << tab << type()
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("constraintDirection") << tab << dir_
       << token::END_STATEMENT << nl;
    os.writeKeyword("translationalAcceleration") << tab << a_
       << token::END_STATEMENT << endl;
}


// ************************************************************************* //
