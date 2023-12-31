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

#include "pairPotential.H"
#include "IOstreams.H"

bool Foam::pairPotential::writeEnergyAndForceTables(Ostream& os) const
{
    Info<< "Writing energy and force tables to file for potential "
        << name_ << endl;

    List< Pair <scalar> > eTab(energyTable());

    List< Pair <scalar> > fTab(forceTable());

    forAll(eTab, e)
    {
        os  << eTab[e].first()
            << token::SPACE
            << eTab[e].second()
            << token::SPACE
            << fTab[e].second()
            << nl;
    }

    return os.good();
}


// ************************************************************************* //
