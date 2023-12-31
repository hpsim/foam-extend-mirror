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

#include "forceCoeffs.H"
#include "objectRegistry.H"
#include "dictionary.H"
#include "foamTime.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forceCoeffs, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forceCoeffs::forceCoeffs
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    forces(name, obr, dict, loadFromFiles),
    liftDir_(vector::zero),
    dragDir_(vector::zero),
    pitchAxis_(vector::zero),
    magUInf_(0.0),
    lRef_(0.0),
    Aref_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::forceCoeffs::~forceCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forceCoeffs::read(const dictionary& dict)
{
    if (active_)
    {
        forces::read(dict);

        // Directions for lift and drag forces, and pitch moment
        dict.lookup("liftDir") >> liftDir_;
        dict.lookup("dragDir") >> dragDir_;
        dict.lookup("pitchAxis") >> pitchAxis_;

        // Free stream velocity magnitude
        dict.lookup("magUInf") >> magUInf_;

        // Reference length and area scales
        dict.lookup("lRef") >> lRef_;
        dict.lookup("Aref") >> Aref_;
    }
}


void Foam::forceCoeffs::writeFileHeader()
{
    if (forcesFilePtr_.valid())
    {
        forcesFilePtr_()
            << "# Time" << tab << "Cd" << tab << "Cl" << tab << "Cm" << endl;
    }
}


void Foam::forceCoeffs::execute()
{
    // Do nothing - only valid on write
}


void Foam::forceCoeffs::end()
{
    // Do nothing - only valid on write
}


void Foam::forceCoeffs::write()
{
    if (active_)
    {
        // Create the forces file if not already created
        makeFile();

        forcesMoments fm = forces::calcForcesMoment();

        scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;

        vector totForce = fm.first().first() + fm.first().second();
        vector totMoment = fm.second().first() + fm.second().second();

        scalar liftForce = totForce & liftDir_;
        scalar dragForce = totForce & dragDir_;
        scalar pitchMoment = totMoment & pitchAxis_;

        scalar Cl = liftForce/(Aref_*pDyn);
        scalar Cd = dragForce/(Aref_*pDyn);
        scalar Cm = pitchMoment/(Aref_*lRef_*pDyn);

        if (Pstream::master())
        {
            forcesFilePtr_()
                << obr_.time().value() << tab
                << Cd << tab << Cl << tab << Cm << endl;

            if (log_)
            {
                Info<< "forceCoeffs output:" << nl
                    << "    Cd = " << Cd << nl
                    << "    Cl = " << Cl << nl
                    << "    Cm = " << Cm << nl
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
