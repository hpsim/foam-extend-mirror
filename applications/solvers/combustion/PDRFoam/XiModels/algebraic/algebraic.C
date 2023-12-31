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

#include "algebraic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(algebraic, 0);
    addToRunTimeSelectionTable(XiModel, algebraic, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::algebraic::algebraic
(
    const dictionary& XiProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su,
    const volScalarField& rho,
    const volScalarField& b,
    const surfaceScalarField& phi
)
:
    XiModel(XiProperties, thermo, turbulence, Su, rho, b, phi),
    XiShapeCoef(readScalar(XiModelCoeffs_.lookup("XiShapeCoef"))),
    XiEqModel_(XiEqModel::New(XiProperties, thermo, turbulence, Su)),
    XiGModel_(XiGModel::New(XiProperties, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiModels::algebraic::~algebraic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::algebraic::Db() const
{
    return XiGModel_->Db();
}


void Foam::XiModels::algebraic::correct()
{
    volScalarField XiEqEta = XiEqModel_->XiEq();
    volScalarField GEta = XiGModel_->G();

    volScalarField R = GEta*XiEqEta/(XiEqEta - 0.999);

    volScalarField XiEqStar = R/(R - GEta);

    Xi_ == 1.0 + (1.0 + (2*XiShapeCoef)*(0.5 - b_))*(XiEqStar - 1.0);
}


bool Foam::XiModels::algebraic::read(const dictionary& XiProperties)
{
    XiModel::read(XiProperties);

    XiModelCoeffs_.lookup("XiShapeCoef") >> XiShapeCoef;

    return true;
}


// ************************************************************************* //
