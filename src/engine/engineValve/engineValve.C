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

#include "engineValve.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::engineValve::adjustCrankAngle(const scalar theta) const
{
    if (theta < liftProfileStart_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta < liftProfileStart_)
        {
            adjustedTheta += liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else if (theta > liftProfileEnd_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta > liftProfileEnd_)
        {
            adjustedTheta -= liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else
    {
        return theta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::engineValve::engineValve
(
    const word& name,
    const polyMesh& mesh,
    const autoPtr<coordinateSystem>& valveCS,
    const word& bottomPatchName,
    const word& poppetPatchName,
    const word& stemPatchName,
    const word& curtainInPortPatchName,
    const word& curtainInCylinderPatchName,
    const word& detachInCylinderPatchName,
    const word& detachInPortPatchName,
    const labelList& detachFaces,
    const graph& liftProfile,
    const scalar minLift,
    const scalar minTopLayer,
    const scalar maxTopLayer,
    const scalar minBottomLayer,
    const scalar maxBottomLayer,
    const scalar diameter
)
:
    name_(name),
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh.time())),
    csPtr_(valveCS),
    bottomPatch_(bottomPatchName, mesh.boundaryMesh()),
    poppetPatch_(poppetPatchName, mesh.boundaryMesh()),
    stemPatch_(stemPatchName, mesh.boundaryMesh()),
    curtainInPortPatch_(curtainInPortPatchName, mesh.boundaryMesh()),
    curtainInCylinderPatch_(curtainInCylinderPatchName, mesh.boundaryMesh()),
    detachInCylinderPatch_(detachInCylinderPatchName, mesh.boundaryMesh()),
    detachInPortPatch_(detachInPortPatchName, mesh.boundaryMesh()),
    detachFaces_(detachFaces),
    liftProfile_(liftProfile),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(minLift),
    minTopLayer_(minTopLayer),
    maxTopLayer_(maxTopLayer),
    minBottomLayer_(minBottomLayer),
    maxBottomLayer_(maxBottomLayer),
    diameter_(diameter)
{}


// Construct from dictionary
Foam::engineValve::engineValve
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh_.time())),
    csPtr_
    (
        coordinateSystem::New
        (
            "coordinateSystem",
            dict.subDict("coordinateSystem")
        )
    ),
    bottomPatch_(dict.lookup("bottomPatch"), mesh.boundaryMesh()),
    poppetPatch_(dict.lookup("poppetPatch"), mesh.boundaryMesh()),
    stemPatch_(dict.lookup("stemPatch"), mesh.boundaryMesh()),
    curtainInPortPatch_
    (
        dict.lookup("curtainInPortPatch"),
        mesh.boundaryMesh()
    ),
    curtainInCylinderPatch_
    (
        dict.lookup("curtainInCylinderPatch"),
        mesh.boundaryMesh()
    ),
    detachInCylinderPatch_
    (
        dict.lookup("detachInCylinderPatch"),
        mesh.boundaryMesh()
    ),
    detachInPortPatch_
    (
        dict.lookup("detachInPortPatch"),
        mesh.boundaryMesh()
    ),
    detachFaces_(dict.lookup("detachFaces")),
    liftProfile_
    (
        "theta",
        "lift",
        name_,
        IFstream
        (
            mesh.time().path()/mesh.time().constant()/
            word(dict.lookup("liftProfileFile"))
        )()
    ),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(readScalar(dict.lookup("minLift"))),
    minTopLayer_(readScalar(dict.lookup("minTopLayer"))),
    maxTopLayer_(readScalar(dict.lookup("maxTopLayer"))),
    minBottomLayer_(readScalar(dict.lookup("minBottomLayer"))),
    maxBottomLayer_(readScalar(dict.lookup("maxBottomLayer"))),
    diameter_(readScalar(dict.lookup("diameter")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::engineValve::lift(const scalar theta) const
{
    return interpolateXY
    (
        adjustCrankAngle(theta),
        liftProfile_.x(),
        liftProfile_.y()
    );
}


bool Foam::engineValve::isOpen() const
{
    return lift(engineDB_.theta()) >= minLift_;
}


Foam::scalar Foam::engineValve::curLift() const
{
    return max
    (
        lift(engineDB_.theta()),
        minLift_
    );
}


Foam::scalar Foam::engineValve::curVelocity() const
{
    return
       -(
             curLift()
           - max
             (
                 lift(engineDB_.theta() - engineDB_.deltaTheta()),
                 minLift_
             )
        )/(engineDB_.deltaT().value() + VSMALL);
}


Foam::labelList Foam::engineValve::movingPatchIDs() const
{
    labelList mpIDs(2);
    label nMpIDs = 0;

    if (bottomPatch_.active())
    {
        mpIDs[nMpIDs] = bottomPatch_.index();
        nMpIDs++;
    }

    if (poppetPatch_.active())
    {
        mpIDs[nMpIDs] = poppetPatch_.index();
        nMpIDs++;
    }

    mpIDs.setSize(nMpIDs);

    return mpIDs;
}


void Foam::engineValve::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK;

    cs().writeDict(os);

    os  << "bottomPatch " << bottomPatch_.name() << token::END_STATEMENT << nl
        << "poppetPatch " << poppetPatch_.name() << token::END_STATEMENT << nl
        << "stemPatch " << stemPatch_.name() << token::END_STATEMENT << nl
        << "curtainInPortPatch " << curtainInPortPatch_.name()
        << token::END_STATEMENT << nl
        << "curtainInCylinderPatch " << curtainInCylinderPatch_.name()
        << token::END_STATEMENT << nl
        << "detachInCylinderPatch " << detachInCylinderPatch_.name()
        << token::END_STATEMENT << nl
        << "detachInPortPatch " << detachInPortPatch_.name()
        << token::END_STATEMENT << nl
        << "detachFaces " << detachFaces_ << token::END_STATEMENT << nl
        << "liftProfile " << nl << token::BEGIN_BLOCK
        << liftProfile_ << token::END_BLOCK << token::END_STATEMENT << nl
        << "minLift " << minLift_ << token::END_STATEMENT << nl
        << "minTopLayer " << minTopLayer_ << token::END_STATEMENT << nl
        << "maxTopLayer " << maxTopLayer_ << token::END_STATEMENT << nl
        << "minBottomLayer " << minBottomLayer_ << token::END_STATEMENT << nl
        << "maxBottomLayer " << maxBottomLayer_ << token::END_STATEMENT << nl
        << "diameter " << diameter_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
