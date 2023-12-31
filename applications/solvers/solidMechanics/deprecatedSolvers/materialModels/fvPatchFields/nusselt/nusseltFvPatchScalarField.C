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

#include "nusseltFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nusseltFvPatchScalarField::nusseltFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    DTName_("undefined"),
    Tinf_(scalar(0)),
    alpha_(p.size(), scalar(0))
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 0;
}


nusseltFvPatchScalarField::nusseltFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    DTName_(dict.lookup("DT")),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    alpha_("alpha", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (Tinf_ < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "unphysical Tinf specified (Tinf = 0 or negative)"
            << exit(FatalError);
    }

    if (min(alpha_) < -SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "unphysical alpha specified (alpha = 0 or negative)" << endl
            << exit(FatalError);
    }
}


nusseltFvPatchScalarField::nusseltFvPatchScalarField
(
    const nusseltFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    DTName_(ptf.DTName_),
    Tinf_(ptf.Tinf_),
    alpha_(ptf.alpha_, mapper)
{}


nusseltFvPatchScalarField::nusseltFvPatchScalarField
(
    const nusseltFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    DTName_(ptf.DTName_),
    Tinf_(ptf.Tinf_),
    alpha_(ptf.alpha_)
{}


nusseltFvPatchScalarField::nusseltFvPatchScalarField
(
    const nusseltFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    DTName_(ptpsf.DTName_),
    Tinf_(ptpsf.Tinf_),
    alpha_(ptpsf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void nusseltFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    alpha_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void nusseltFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const nusseltFvPatchScalarField& npsf =
        refCast<const nusseltFvPatchScalarField>(ptf);

    alpha_.rmap(npsf.alpha_, addr);
}


// Update the coefficients associated with the patch field
void nusseltFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField Tinternal = patchInternalField();

    // Lookup temperature diffusivity of the patch
    const fvPatchField<scalar>& DT =
        this->patch().lookupPatchField<volScalarField, scalar>(DTName_);

    // Calculate flux
    scalarField tempFlux = alpha_*(Tinternal - Tinf_);

    refValue() =
        neg(tempFlux)*
        min
        (
            Tinternal - tempFlux/(DT*patch().deltaCoeffs()),
            Tinf_
        )
      + pos(tempFlux)*
        max
        (
            Tinternal - tempFlux/(DT*patch().deltaCoeffs()),
            Tinf_
        );

    refGrad() = -tempFlux;
    valueFraction() = pos(tempFlux);

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void nusseltFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("DT") << DTName_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << endl;
    alpha_.writeEntry("alpha", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nusseltFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
