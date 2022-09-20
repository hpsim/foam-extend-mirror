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

Class
    fluidThermalModel

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fluidThermalModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fluidThermalModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermalModel::fluidThermalModel
(
    const dictionary& dict,
    const volScalarField& T
)
:
    dict_(dict),
    T_(T)
{
    // If thermal dictionary is found, read it
    if (dict_.found("thermal"))
    {
        // Get thermal law
        lawPtr_ = thermalLaw::New("law", T_, dict_.subDict("thermal"));

        PtrList<entry> entries(dict_.subDict("thermal").lookup("sources"));
        sources_.setSize(entries.size());

        forAll (sources_, sourceI)
        {
            sources_.set
            (
                sourceI,
                thermalSource::New
                (
                    entries[sourceI].keyword(),
                    T,
                    entries[sourceI].dict()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermalModel::S() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "heatSource",
                T_.time().timeName(),
                T_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionedScalar
            (
                "zero",
                // Watt/m^3 divided by rho*Cp
                dimEnergy/dimTime/dimVolume,
                scalar(0)
            )
        )
    );
    volScalarField& source = tsource();

    forAll(sources_, sourceI)
    {
        sources_[sourceI].addSource(source);
    }

    return tsource;
}


// ************************************************************************* //
