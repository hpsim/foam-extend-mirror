/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
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

#include "triSurfaceDetectFeatureEdges.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceDetectFeatureEdges::triSurfaceDetectFeatureEdges
(
    triSurf& surface,
    const scalar angleDeviation
)
:
    surf_(surface),
    featureEdges_(surf_.edges().size(), direction(0)),
    angleTolerance_(angleDeviation)
{
    if( Pstream::parRun() )
        FatalError << "Feature edges detection does not run in parallel"
            << exit(FatalError);

    detectFeatureEdgesAngleCriterion();
}

triSurfaceDetectFeatureEdges::~triSurfaceDetectFeatureEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectFeatureEdges::detectFeatureEdges()
{
    const edgeLongList& edges = surf_.edges();
    triSurfModifier surfMod(surf_);
    edgeLongList& featureEdges = surfMod.featureEdgesAccess();
    featureEdges.clear();

    forAll(featureEdges_, eI)
    {
        if( featureEdges_[eI] )
            featureEdges.append(edges[eI]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
