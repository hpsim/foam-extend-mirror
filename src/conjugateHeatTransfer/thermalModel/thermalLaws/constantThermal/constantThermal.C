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

#include "constantThermal.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantThermal, 0);
    addToRunTimeSelectionTable(thermalLaw, constantThermal, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::constantThermal::constantThermal
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    thermalLaw(name, T, dict),
    rho_(dict.lookup("rho")),
    C_(dict.lookup("C")),
    k_(dict.lookup("k")),
    alpha_(dict.lookup("alpha")),
    T0_(dict.lookup("T0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantThermal::~constantThermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constantThermal::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::constantThermal::C() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "C",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            C_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::constantThermal::k() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "kTmp",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            k_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::constantThermal::alpha() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "alpha",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            alpha_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::constantThermal::T0() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "T0",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            T0_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


// ************************************************************************* //
