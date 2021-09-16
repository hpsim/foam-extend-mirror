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

Description

\*---------------------------------------------------------------------------*/

#include "wallFvPatch.H"
#include "fvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(wallFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, wallFvPatch, polyPatch);

} // End namespace Foam


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::wallFvPatch::makeCorrVecs(fvsPatchVectorField& cv) const
{
    cv = vector::zero;
}


void Foam::wallFvPatch::updatePhi
(
    DimensionedField<scalar, volMesh>& V,
    DimensionedField<scalar, volMesh>& V0,
    surfaceScalarField& phi
) const
{
    if (wallPolyPatch_.closedSolidBodyMotion())
    {
        const scalar phiAvg = gAverage(phi.boundaryField()[index()]);

        phi.boundaryField()[index()] -= phiAvg;

        Info<< "wallFvPatch::updatePhi updating average for closed wall: "
            << phiAvg << " corrected: "
            << gAverage(phi.boundaryField()[index()])
            << endl;
    }
}


// ************************************************************************* //
