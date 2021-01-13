/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved
    Refactoring by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#include "chtRcThermalDiffusivityResistanceFvPatchScalarField.H"
#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "harmonic.H"
#include "radiationConstants.H"
#include "VectorN.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(p, iF),
    conductivity_(p.size(), 0)
{}


Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRegionCoupleBase(p, iF, dict),
    conductivity_("conductivity", dict, p.size())
{}


Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const chtRcThermalDiffusivityResistanceFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRegionCoupleBase(ptf, p, iF, mapper),
    conductivity_(ptf.conductivity_)
{}


Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::
chtRcThermalDiffusivityResistanceFvPatchScalarField
(
    const chtRcThermalDiffusivityResistanceFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(ptf, iF),
    conductivity_(ptf.conductivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    calcThermalDiffusivity(*this, shadowPatchField());
}


void
Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::
calcThermalDiffusivity
(
    chtRegionCoupleBase& owner,
    const chtRegionCoupleBase& neighbour
) const
{
    const fvPatch& p = owner.patch();
    const fvMesh& mesh = p.boundaryMesh().mesh();

    const chtRcTemperatureFvPatchScalarField& TwOwn =
        dynamic_cast<const chtRcTemperatureFvPatchScalarField&>
        (
            p.lookupPatchField<volScalarField, scalar>("T")
        );
    scalarField& k = owner;
    const scalarField& fOwn = owner.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());
    scalarField TwNei(p.size());

    scalarField QrOwn(p.size(), 0.0);
    scalarField fourQroOwn(p.size(), 0.0);
    scalarField fourQroNei(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        QrOwn += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQroOwn += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    scalarField Qr = QrOwn;
    scalarField cond(p.size());

    {
        Field<VectorN<scalar, 5> > lData(neighbour.size());

        const scalarField& lfNei = neighbour.originalPatchField();
        scalarField lTcNei = TwOwn.shadowPatchField().patchInternalField();
        const scalarField& lTwNei = TwOwn.shadowPatchField();

        forAll(lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
            lData[facei][2] = lTwNei[facei];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            const scalarField& lQrNei =
                owner.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll(lData, facei)
            {
                lData[facei][3] = lQrNei[facei];
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(owner))
        {
            forAll(lData, facei)
            {
                lData[facei][4] = conductivity_[facei];
            }
        }

        const Field<VectorN<scalar, 5> > iData =
            owner.regionCouplePatch().interpolate(lData);

        forAll(iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
            TwNei[facei] = iData[facei][2];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            forAll(iData, facei)
            {
                fourQroNei[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
                Qr[facei] += iData[facei][3];
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(owner))
        {
            forAll(lData, facei)
            {
                cond[facei] = iData[facei][4];
            }
        }
        else
        {
            cond = conductivity_;
        }
    }

    // Do interpolation
    harmonic<scalar> interp(mesh);
    const scalarField weights = interp.weights(fOwn, fNei, p);
    scalarField kHarm = weights*fOwn + (1.0 - weights)*fNei;
    kHarm *= cond/(kHarm*p.deltaCoeffs() + cond);

    const scalarField kOwn = fOwn/((1.0 - p.weights())*p.magLongDeltas());
    const scalarField kNei = fNei/(p.weights()*p.magLongDeltas());

    scalarField temp = TwNei*(cond + kNei) + fourQroNei;

    // Expression is equivalent to the one above
    k = kOwn*
    (
        TwOwn*
        (
            TwNei*
            (
                cond*(kNei*(TcOwn - TcNei) - fourQroOwn - fourQroNei - Qr)
              - kNei*(fourQroOwn + QrOwn)
            )
          + fourQroNei*(TcOwn*cond - fourQroOwn - QrOwn)
        )
      + (TcOwn*fourQroOwn*temp)
    );

    k /= stabilise
    (
        (TcOwn - TcNei)*
        (
            TwOwn*
            (
                kOwn*temp + cond*(TwNei*kNei + fourQroNei)
            )
          + fourQroOwn*temp
        ),
        SMALL
    );

    k /= p.deltaCoeffs();

    forAll(k, facei)
    {
        k[facei] = max(min(k[facei], 100*kHarm[facei]), 0.01*kHarm[facei]);
    }

    owner.fvPatchScalarField::updateCoeffs();
}


void
Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::calcTemperature
(
    chtRcTemperatureFvPatchScalarField& TwOwn,
    const chtRcTemperatureFvPatchScalarField& neighbour,
    const chtRegionCoupleBase& ownerK
) const
{
    const fvPatch& p = TwOwn.patch();

    const scalarField& fOwn = ownerK.originalPatchField();
    const scalarField TcOwn = TwOwn.patchInternalField();

    scalarField fNei(p.size());
    scalarField TcNei(p.size());
    scalarField TwNei(p.size());

    scalarField QrOwn(p.size(), 0.0);
    scalarField fourQroOwn(p.size(), 0.0);
    scalarField fourQroNei(p.size(), 0.0);

    if (TwOwn.radiation())
    {
        QrOwn += p.lookupPatchField<volScalarField, scalar>("Qr");
        fourQroOwn += 4.0*radiation::sigmaSB.value()*pow4(TwOwn);
    }

    scalarField Qr = QrOwn;
    scalarField cond(p.size());

    {
        Field<VectorN<scalar, 5> > lData(neighbour.size());

        const scalarField& lfNei =
            ownerK.shadowPatchField().originalPatchField();
        scalarField lTcNei =
            TwOwn.shadowPatchField().patchInternalField();
        const scalarField& lTwNei =
            TwOwn.shadowPatchField();

        forAll(lData, facei)
        {
            lData[facei][0] = lTcNei[facei];
            lData[facei][1] = lfNei[facei];
            lData[facei][2] = lTwNei[facei];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            const scalarField& lQrNei =
                TwOwn.lookupShadowPatchField<volScalarField, scalar>("Qr");

            forAll(lData, facei)
            {
                lData[facei][3] = lQrNei[facei];
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(ownerK))
        {
            forAll(lData, facei)
            {
                lData[facei][4] = conductivity_[facei];
            }
        }

        const Field<VectorN<scalar, 5> > iData =
            TwOwn.regionCouplePatch().interpolate(lData);

        forAll(iData, facei)
        {
            TcNei[facei] = iData[facei][0];
            fNei[facei] = iData[facei][1];
            TwNei[facei] = iData[facei][2];
        }

        if(TwOwn.shadowPatchField().radiation())
        {
            forAll(iData, facei)
            {
                Qr[facei] += iData[facei][3];
                fourQroNei[facei] +=
                    4.0*radiation::sigmaSB.value()*pow4(iData[facei][2]);
            }
        }

        if(isA<chtRcThermalDiffusivitySlaveFvPatchScalarField>(ownerK))
        {
            forAll(lData, facei)
            {
                cond[facei] = iData[facei][4];
            }
        }
        else
        {
            cond = conductivity_;
        }
    }

    const scalarField kOwn = fOwn/((1.0 - p.weights())*p.magLongDeltas());
    const scalarField kNei = fNei/(p.weights()*p.magLongDeltas());

    scalarField temp = fourQroOwn + QrOwn + kOwn*TcOwn;

    TwOwn *=
    (
        TwNei*
        (
            cond*
            (
                kNei*TcNei + fourQroOwn + fourQroNei + Qr + kOwn*TcOwn
            )
          + kNei*temp
        )
      + fourQroNei*temp
    )
    /
    (
        TwOwn*
        (
            TwNei*(cond*(kOwn + kNei) + kOwn*kNei)
          + fourQroNei*(kOwn + cond)
        )
      + fourQroOwn*(TwNei*(cond + kNei) + fourQroNei)
    );

    TwOwn.fvPatchScalarField::updateCoeffs();
}


void Foam::chtRcThermalDiffusivityResistanceFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    conductivity_.writeEntry("conductivity", os);
    os.writeKeyword("remoteField")
        << remoteFieldName() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        chtRcThermalDiffusivityResistanceFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
