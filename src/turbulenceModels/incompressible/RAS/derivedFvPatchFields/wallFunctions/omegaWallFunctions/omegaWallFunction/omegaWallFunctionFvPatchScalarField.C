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

#include "omegaWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void omegaWallFunctionFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("omegaWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("RASModel::G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075)
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_)
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075))
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf, iF),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(GName_))
    {
        InfoIn("void omegaWallFunctionFvPatchScalarField::updateCoeffs()")
            << "Cannot access " << GName_ << " field for patch "
            << patch().name() << ".  Evaluating as zeroGradient"
            << endl;

        fvPatchScalarField::updateCoeffs();
        zeroGradientFvPatchScalarField::evaluate();

        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    // Note: omega is now a refValue and set in fixedInternalValueFvPatchField
    // HJ, 3/Aug/2011
    scalarField& omega = refValue();

    const scalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& nuw =
        lookupPatchField<volScalarField, scalar>(nuName_);

    const scalarField& nutw =
        lookupPatchField<volScalarField, scalar>(nutName_);

    const fvPatchVectorField& Uw =
        lookupPatchField<volVectorField, vector>(UName_);

    vectorField n = patch().nf();

    const scalarField magGradUw = mag(Uw.snGrad());

    const labelList& faceCells = patch().faceCells();

    // Set omega and G
    forAll(nutw, faceI)
    {
        label faceCellI = faceCells[faceI];

        scalar yPlus = Cmu25*y[faceI]*sqrt(k[faceCellI])/nuw[faceI];

        scalar omegaVis = 6.0*nuw[faceI]/(beta1_*sqr(y[faceI]));

        scalar omegaLog = sqrt(k[faceCellI])/(Cmu25*kappa_*y[faceI]);

        // Note: omega is now a refValue and set in
        // fixedInternalValueFvPatchField
        // HJ, 3/Aug/2011
        omega[faceI] = sqrt(sqr(omegaVis) + sqr(omegaLog));

        if (yPlus > yPlusLam)
        {
            G[faceCellI] =
                (nutw[faceI] + nuw[faceI])
               *magGradUw[faceI]
               *Cmu25*sqrt(k[faceCellI])
               /(kappa_*y[faceI]);
        }
        else
        {
            G[faceCellI] = 0.0;
        }
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}


void omegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
