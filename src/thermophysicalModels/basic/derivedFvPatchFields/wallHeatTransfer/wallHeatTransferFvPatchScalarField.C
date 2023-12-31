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

#include "wallHeatTransferFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatTransferFvPatchScalarField::wallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_(p.size(), 0.0),
    alphaWall_(p.size(), 0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::wallHeatTransferFvPatchScalarField::wallHeatTransferFvPatchScalarField
(
    const wallHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Tinf_(ptf.Tinf_, mapper),
    alphaWall_(ptf.alphaWall_, mapper)
{}


Foam::wallHeatTransferFvPatchScalarField::wallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_("Tinf", dict, p.size()),
    alphaWall_("alphaWall", dict, p.size())
{
    refValue() = Tinf_;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::wallHeatTransferFvPatchScalarField::wallHeatTransferFvPatchScalarField
(
    const wallHeatTransferFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    Tinf_(tppsf.Tinf_),
    alphaWall_(tppsf.alphaWall_)
{}


Foam::wallHeatTransferFvPatchScalarField::wallHeatTransferFvPatchScalarField
(
    const wallHeatTransferFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    Tinf_(tppsf.Tinf_),
    alphaWall_(tppsf.alphaWall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallHeatTransferFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    Tinf_.autoMap(m);
    alphaWall_.autoMap(m);
}


void Foam::wallHeatTransferFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const wallHeatTransferFvPatchScalarField& tiptf =
        refCast<const wallHeatTransferFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);
    alphaWall_.rmap(tiptf.alphaWall_, addr);
}


void Foam::wallHeatTransferFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );

    const label patchi = patch().index();

    const scalarField& Tw = thermo.T().boundaryField()[patchi];
    scalarField Cpw = thermo.Cp(Tw, patchi);

    valueFraction() =
        1.0/
        (
            1.0
          + Cpw*thermo.alpha().boundaryField()[patchi]
           *patch().deltaCoeffs()/alphaWall_
        );

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::wallHeatTransferFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    Tinf_.writeEntry("Tinf", os);
    alphaWall_.writeEntry("alphaWall", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, wallHeatTransferFvPatchScalarField);
}

// ************************************************************************* //
