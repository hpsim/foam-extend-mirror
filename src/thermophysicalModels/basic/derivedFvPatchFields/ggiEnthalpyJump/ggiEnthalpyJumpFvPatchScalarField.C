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

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Gregor Cvijetic, FMENA Zagreb.
    Hrvoje Jasak, Wikki Ltd.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

\*---------------------------------------------------------------------------*/

#include "ggiEnthalpyJumpFvPatchScalarField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ggiEnthalpyJumpFvPatchScalarField::ggiEnthalpyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpGgiFvPatchScalarField(p, iF),
    rotating_(false),
    jump_(this->size(), pTraits<scalar>::zero)
{
    // Take rotating_ from temperature BC
    if (iF.name() == db().mangleFileName("h"))
    {

        const fvPatchScalarField& T =
            this->lookupPatchField
                (
                    "T",
                    reinterpret_cast<const volScalarField*>(0),
                    reinterpret_cast<const scalar*>(0)
                );

        const ggiEnthalpyJumpFvPatchScalarField& Tpatch =
            dynamic_cast<const ggiEnthalpyJumpFvPatchScalarField&>(T);

        rotating_ = Tpatch.rotating();
    }
}


ggiEnthalpyJumpFvPatchScalarField::ggiEnthalpyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    jumpGgiFvPatchScalarField(p, iF),
    rotating_(dict.lookup("rotating")),
    jump_(this->size(), pTraits<scalar>::zero)
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::blocking);
    }
}


ggiEnthalpyJumpFvPatchScalarField::ggiEnthalpyJumpFvPatchScalarField
(
    const ggiEnthalpyJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpGgiFvPatchScalarField(ptf, p, iF, mapper),
    rotating_(ptf.rotating_),
    jump_(ptf.jump_, mapper)
{}


ggiEnthalpyJumpFvPatchScalarField::ggiEnthalpyJumpFvPatchScalarField
(
    const ggiEnthalpyJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    jumpGgiFvPatchScalarField(ptf, iF),
    rotating_(ptf.rotating_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ggiEnthalpyJumpFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpGgiFvPatchScalarField::autoMap(m);
    jump_.autoMap(m);
}


void ggiEnthalpyJumpFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    jumpGgiFvPatchScalarField::rmap(ptf, addr);

    // rmap jump
    const ggiEnthalpyJumpFvPatchScalarField& ejPtf =
        refCast<const ggiEnthalpyJumpFvPatchScalarField>(ptf);

    jump_.rmap(ejPtf.jump_, addr);
}


void Foam::ggiEnthalpyJumpFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get access to relative and rotational velocity
    const word URotName("URot");
    const word UThetaName("UTheta");

    jump_ = 0;

    // If this is rothalpy, calculate the jump, otherwise leave it equal to 0.
    if
    (
        dimensionedInternalField().name() == db().mangleFileName("i")
    )
    {
        if
        (
            !this->db().objectRegistry::found(URotName)
         || !this->db().objectRegistry::found(UThetaName)
        )
        {
            // Velocities not available, do not update
            InfoIn
            (
                "void gradientEnthalpyFvPatchScalarField::"
                "updateCoeffs(const vectorField& Up)"
            )   << "Velocity fields " << URotName << " or "
                << UThetaName << " not found.  "
                << "Performing enthalpy value update" << endl;

            jump_ = 0;
        }
        else
        {
            const fvPatchVectorField& URotp =
                lookupPatchField<volVectorField, vector>(URotName);

            const fvPatchScalarField& UThetap =
                lookupPatchField<volScalarField, scalar>(UThetaName);

            if (rotating_)
            {
                // We can either make jump_ on neighbour field and interpolate
                // (in jumpOverlapGgi) or interpolate first, then add jump_ for
                // internalField
                const scalarField UThetaIn = UThetap.patchInternalField();

                jump_ = -mag(UThetaIn)*mag(URotp.patchInternalField());
            }
            else
            {
                // Patch neighbour field neccessary as URot = 0 on stator
                const scalarField UThetaIn = UThetap.patchNeighbourField();

                jump_ = mag(UThetaIn)*mag(URotp.patchNeighbourField());
            }
        }
    }

    jumpGgiFvPatchScalarField::updateCoeffs();
}


void ggiEnthalpyJumpFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("patchType")
        << ggiFvPatch::typeName << token::END_STATEMENT << nl;
    os.writeKeyword("rotating")
        << rotating_ << token::END_STATEMENT << nl;

    IOstream::streamFormat fmt0 = os.format(IOstream::ASCII);
    os.format(fmt0);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    ggiEnthalpyJumpFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
