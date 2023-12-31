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
    verticalValvesGambit

\*---------------------------------------------------------------------------*/

#include "pistonSliding.H"
#include "layerAdditionRemoval.H"
#include "attachDetach.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "engineTime.H"
#include "pointField.H"
#include "fvPatchField.H"
#include "Switch.H"
#include "symmetryFvPatch.H"
#include "tetMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pistonSliding, 0);
    addToRunTimeSelectionTable(engineTopoChangerMesh, pistonSliding, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




bool Foam::pistonSliding::realDeformation() const
{

    bool deformationValve = false;
    forAll(valves(), valveI)
    {

        scalar maxLayer = piston().minLayer();

        if(valves()[valveI].bottomPatchID().active())
        {
            maxLayer = max(maxLayer, valves()[valveI].minBottomLayer());
        }

        scalar valveDisplacement = valves_[valveI].curVelocity()*valves_[valveI].cs().axis().z()*engTime().deltaT().value()  ;
        if(valvePistonPosition()[valveI] + engTime().pistonDisplacement().value() >
        valveBottomPosition_[valveI] + valveDisplacement - 5.0*maxLayer - 0.001 )
        {
            deformationValve = true;
        }
    }

    if(deformationValve)
    {
        return true;
    }
    else  if (virtualPistonPosition() + engTime().pistonDisplacement().value() > deckHeight_ - SMALL)
    {
        return true;
    }
    else
    {
        return deformation();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pistonSliding::pistonSliding
(
    const IOobject& io
)
:
    engineTopoChangerMesh(io),
    piston_(*this, engTime().engineDict().subDict("piston")),
    valves_(*this, engTime().engineDict().lookup("pistonSliding")),
    movingPointsMaskTopPtr_(nullptr),
    movingPointsMaskBottomPtr_(nullptr),
    deformSwitch_(readScalar(engTime().engineDict().lookup("deformAngle"))),
    valveTopTol_(readScalar(engTime().engineDict().lookup("valveTopTol"))),
    pistonPosition_(-GREAT),
    virtualPistonPosition_(-GREAT),
    valveTopPosition_(nValves(),-GREAT),
    valveBottomPosition_(nValves(),GREAT),
    valvePistonPosition_(nValves(),GREAT),
    deckHeight_(GREAT),
    minValveZ_(nValves()),
    poppetValveTol_(readScalar(engTime().engineDict().lookup("poppetValveTol"))),
    bottomValveTol_(readScalar(engTime().engineDict().lookup("bottomValveTol"))),
    msPtr_(motionSolver::New(*this)),
    isReallyClosed_(valves().size(), false),
    correctPointsMotion_(engTime().engineDict().lookup("correctPointsMotion"))


{
    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pistonSliding::setBoundaryVelocity(volVectorField& U)
{


    // Set valve velociaty
    forAll (valves(), valveI)
    {

        vector valveVel =
            valves()[valveI].curVelocity()*valves()[valveI].cs().axis();

        // If valve is present in geometry, set the motion
        if (valves()[valveI].curtainInPortPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].curtainInPortPatchID().index()] ==
//                valveVel;
                vector::zero;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].curtainInCylinderPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].curtainInCylinderPatchID().index()] ==
//                valveVel;
                vector::zero;
        }

        // If valve is present in geometry, set the motion
        if (valves()[valveI].poppetPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].poppetPatchID().index()] ==
                valveVel;
        }

        if (valves()[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].bottomPatchID().index()] ==
                valveVel;
        }


        // If valve is present in geometry, set the motion
        if (valves()[valveI].stemPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            U.boundaryField()[valves()[valveI].stemPatchID().index()] ==
                valveVel;
        }


    }

}

bool Foam::pistonSliding::inValve(const point& p, const label& i) const
{
    scalar valveX = valves_[i].cs().origin().x();
    scalar valveY = valves_[i].cs().origin().y();
    return (sqrt(sqr(p.x()-valveX)+sqr(p.y()-valveY)) < 0.5*valves_[i].diameter());
}

bool Foam::pistonSliding::inPiston(const point& p) const
{
    scalar pistonX = piston_.cs().origin().x();
    scalar pistonY = piston_.cs().origin().y();
    return (sqrt(sqr(p.x()-pistonX)+sqr(p.y()-pistonY)) < 0.5*engTime().bore().value());
}


bool Foam::pistonSliding::isACylinderHeadFace
(
    const labelList& cylHeadFaces,
    const label face
)
{
    forAll(cylHeadFaces, i)
    {
        if(cylHeadFaces[i] == face)
        {
            return true;
        }
    }

    return false;
}



// ************************************************************************* //
