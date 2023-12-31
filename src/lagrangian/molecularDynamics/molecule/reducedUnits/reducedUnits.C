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

#include "reducedUnits.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const  Foam::scalar Foam::reducedUnits::kb = 1.3806504e-23;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reducedUnits::calcRefValues()
{
    if
    (
        refTime_ < VSMALL
     || refLength_ < VSMALL
     || refMass_ < VSMALL
    )
    {
        FatalErrorIn("Foam::reducedUnits::calcRefValues() ")
            << "One of more referencence values too small for floating point "
            << "calculation: "
            << "refTime_ = " << refTime_
            << ", refLength = " << refTemp_
            << ", refMass = " << refMass_
            << nl << abort(FatalError);
    }

    refEnergy_ = refLength_*refLength_*refMass_/(refTime_*refTime_);

    refTemp_ = refEnergy_ / kb;

    refForce_ = refEnergy_/refLength_;

    refVelocity_ = Foam::sqrt(refEnergy_/refMass_);

    refVolume_ = Foam::pow(refLength_,3.0);

    refPressure_ = refEnergy_/refVolume_;

    refMassDensity_ = refMass_/refVolume_;

    refNumberDensity_ = 1.0/refVolume_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reducedUnits::reducedUnits()
:
    refLength_(1e-9),
    refTime_(1e-12),
    refMass_(1.660538782e-27)
{
    calcRefValues();
}


Foam::reducedUnits::reducedUnits
(
    scalar refLength,
    scalar refTime,
    scalar refMass
)
:
    refLength_(refLength),
    refTime_(refTime),
    refMass_(refMass)
{
    calcRefValues();
}


Foam::reducedUnits::reducedUnits(const IOdictionary& reducedUnitsDict)
:
    refLength_(),
    refTime_(),
    refMass_()
{
    setRefValues(reducedUnitsDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reducedUnits::~reducedUnits()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reducedUnits::setRefValues
(
    scalar refLength,
    scalar refTime,
    scalar refMass
)
{
    refLength_ = refLength;

    refTime_ = refTime;

    refMass_ = refMass;

    calcRefValues();
}


void Foam::reducedUnits::setRefValues
(
    const IOdictionary& reducedUnitsDict
)
{
    refLength_ = readScalar(reducedUnitsDict.lookup("refLength"));

    refTime_ = readScalar(reducedUnitsDict.lookup("refTime"));

    refMass_  = readScalar(reducedUnitsDict.lookup("refMass"));

    calcRefValues();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::reducedUnits::operator=(const reducedUnits& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::reducedUnits::operator=(const Foam::reducedUnits&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// ************************************************************************* //
