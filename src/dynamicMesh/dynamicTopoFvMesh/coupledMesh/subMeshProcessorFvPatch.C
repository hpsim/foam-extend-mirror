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

Class
    subMeshProcessorFvPatch

Description
    Functions specific to subMeshProcessorFvPatch

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "subMeshProcessorFvPatch.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(subMeshProcessorFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, subMeshProcessorFvPatch, polyPatch);


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void subMeshProcessorFvPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


void subMeshProcessorFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    dc = 1.0/(nf() & fvPatch::delta());
}


tmp<vectorField> subMeshProcessorFvPatch::delta() const
{
    return fvPatch::delta();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
