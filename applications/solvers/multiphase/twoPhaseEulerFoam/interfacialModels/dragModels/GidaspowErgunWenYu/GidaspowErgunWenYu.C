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

#include "GidaspowErgunWenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        GidaspowErgunWenYu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GidaspowErgunWenYu::GidaspowErgunWenYu
(
    const dictionary& interfaceDict,
    const volScalarField& alpha,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    dragModel(interfaceDict, alpha, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GidaspowErgunWenYu::K
(
    const volScalarField& Ur
) const
{
    volScalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));

    volScalarField bp = pow(beta, -2.65);
    volScalarField Re = max(Ur*phasea_.d()/phaseb_.nu(), scalar(1.0e-3));

    volScalarField Cds = 24.0*(1.0 + 0.15*pow(Re, 0.687))/Re;

    forAll(Re, celli)
    {
        if(Re[celli] > 1000.0)
        {
            Cds[celli] = 0.44;
        }
    }

    // Wen and Yu (1966)
    tmp<volScalarField> tKWenYu = 0.75*Cds*phaseb_.rho()*Ur*bp/phasea_.d();
    volScalarField& KWenYu = tKWenYu();

    // Ergun
    forAll (beta, cellj)
    {
        if (beta[cellj] <= 0.8)
        {
            KWenYu[cellj] =
                150.0*alpha_[cellj]*phaseb_.nu().value()*phaseb_.rho().value()
               /sqr(beta[cellj]*phasea_.d().value())
              + 1.75*phaseb_.rho().value()*Ur[cellj]
               /(beta[cellj]*phasea_.d().value());
        }
    }

    return tKWenYu;
}


// ************************************************************************* //
