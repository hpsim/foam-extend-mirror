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

#include "basicThermo.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnthalpyFvPatchScalarField.H"
#include "gradientEnthalpyFvPatchScalarField.H"
#include "mixedEnthalpyFvPatchScalarField.H"
#include "fixedInternalEnergyFvPatchScalarField.H"
#include "gradientInternalEnergyFvPatchScalarField.H"
#include "mixedInternalEnergyFvPatchScalarField.H"

/* * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicThermo::hBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnthalpyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnthalpyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnthalpyFvPatchScalarField::typeName;
        }
    }

    return hbt;
}


void Foam::basicThermo::hBoundaryCorrection(volScalarField& h)
{
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnthalpyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnthalpyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


Foam::wordList Foam::basicThermo::eBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList ebt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            ebt[patchi] = fixedInternalEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            ebt[patchi] = gradientInternalEnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            ebt[patchi] = mixedInternalEnergyFvPatchScalarField::typeName;
        }
    }

    return ebt;
}


void Foam::basicThermo::eBoundaryCorrection(volScalarField& e)
{
    volScalarField::GeometricBoundaryField& ebf = e.boundaryField();

    forAll(ebf, patchi)
    {
        if (isA<gradientInternalEnergyFvPatchScalarField>(ebf[patchi]))
        {
            refCast<gradientInternalEnergyFvPatchScalarField>(ebf[patchi])
                .gradient() = ebf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedInternalEnergyFvPatchScalarField>(ebf[patchi]))
        {
            refCast<mixedInternalEnergyFvPatchScalarField>(ebf[patchi])
                .refGrad() = ebf[patchi].fvPatchField::snGrad();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicThermo::basicThermo(const fvMesh& mesh, const objectRegistry& obj)
:
    IOdictionary
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            obj,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            obj,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    psi_
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0, 0, 0)
    ),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            obj,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    alpha_
    (
        IOobject
        (
            "alpha",
            mesh.time().timeName(),
            obj,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::psi() const
{
    return psi_;
}


Foam::volScalarField& Foam::basicThermo::h()
{
    notImplemented("basicThermo::h()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::basicThermo::h() const
{
    notImplemented("basicThermo::h() const");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "basicThermo::h"
        "(const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::h
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "basicThermo::h"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::volScalarField& Foam::basicThermo::hs()
{
    notImplemented("basicThermo::hs()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::basicThermo::hs() const
{
    notImplemented("basicThermo::hs() const");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "basicThermo::hs"
        "(const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "basicThermo::hs"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::basicThermo::hc() const
{
    notImplemented("basicThermo::hc()");
    return volScalarField::null();
}


Foam::volScalarField& Foam::basicThermo::e()
{
    notImplemented("basicThermo::e()");
    return const_cast<volScalarField&>(volScalarField::null());
}


const Foam::volScalarField& Foam::basicThermo::e() const
{
    notImplemented("basicThermo::e()");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::e
(
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "basicThermo::e"
        "(const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::e
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "basicThermo::e"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(nullptr);
}


const Foam::volScalarField& Foam::basicThermo::T() const
{
    return T_;
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "basicThermo::Cp"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::Cp
(
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "basicThermo::Cp"
        "(const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::basicThermo::Cp() const
{
    notImplemented("basicThermo::Cp() const");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "basicThermo::Cv"
        "(const scalarField& T, const label patchi) const"
    );
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::basicThermo::Cv() const
{
    notImplemented("basicThermo::Cv() const");
    return volScalarField::null();
}


const Foam::volScalarField& Foam::basicThermo::mu() const
{
    return mu_;
}


const Foam::volScalarField& Foam::basicThermo::alpha() const
{
    return alpha_;
}


bool Foam::basicThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
