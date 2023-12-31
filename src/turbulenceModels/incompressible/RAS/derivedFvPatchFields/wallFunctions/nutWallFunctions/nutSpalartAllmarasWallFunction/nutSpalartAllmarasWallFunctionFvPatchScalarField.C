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

#include "nutSpalartAllmarasWallFunctionFvPatchScalarField.H"
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
nutSpalartAllmarasWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField magGradU = mag(Uw.snGrad());
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    return Foam::max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


tmp<scalarField> nutSpalartAllmarasWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];

    const fvPatchVectorField& Uw =
        turbModel.U().boundaryField()[patchI];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);

    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];
    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();

    forAll(uTau, faceI)
    {
        const scalar& magUpara = magUp[faceI];

        scalar ut = sqrt((nutw[faceI] + nuw[faceI])*magGradU[faceI]);

        if (ut > VSMALL)
        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                const scalar kUu = min(kappa_*magUpara/ut, 50);
                const scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                const scalar f =
                    - ut*y[faceI]/nuw[faceI]
                    + magUpara/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                const scalar df =
                    y[faceI]/nuw[faceI]
                  + magUpara/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                const scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > VSMALL && err > 0.01 && ++iter < 10);

            uTau[faceI] = max(0.0, ut);
        }
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
nutSpalartAllmarasWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutSpalartAllmarasWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
