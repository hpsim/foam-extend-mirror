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
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "attachDetach.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pistonSliding::addZonesAndModifiers()
{
    // Add the zones and mesh modifiers to operate piston motion

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
     || topoChanger_.size() > 0
    )
    {
        Info<< "Time = " << engTime().theta() << endl;
        Info<< "void Foam::verticalValvesGambit::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        setVirtualPositions();
//        checkAndCalculate();

        Info << "Point zones found = " << pointZones().size() << endl;
        Info << "Face zones found = " << faceZones().size() << endl;
        Info << "Cell zones found = " << cellZones().size() << endl;

        return;

    }

    if
    (
        engTime().engineDict().found("zOffsetGambit")
     && engTime().engineDict().found("zDisplGambit")
    )
    {
        Info << "Assembling the cylinder mesh" << endl;

        scalar zOffset
        (
            readScalar(engTime().engineDict().lookup("zOffsetGambit"))
        );

        scalar zDispl
        (
            readScalar(engTime().engineDict().lookup("zDisplGambit"))
        );

        pointField pDispl = points();

        forAll(points(), pointI)
        {
            const point p = points()[pointI];

            if(p.z() >= zOffset)
            {
                pDispl[pointI].z() -= zDispl;
            }
        }

        movePoints(pDispl);
        write();
        resetMotion();

        Info << "Cylinder mesh assembled" << endl;
    }

    Info << "checkAndCalculate()" << endl;
//    checkAndCalculate();

    Info<< "Time = " << engTime().theta() << endl
        << "Adding zones to the engine mesh" << endl;

/*
    Point zones
    1) Piston points
    1) Cut point zone for liner in head

    nValves*
    1) cutPointsV
    2) valveTopPoints
    3) valveBottomPoints
*/
    DynamicList<pointZone*> pz;

/*
    Face zones
    1) Piston layer faces

    nValves*
    1) valveTopLayerFaces
    2) valveBottomLayerFaces
    3) valveCurtainPort
    4) valveCurtainCyl
    5) cutFaceV
*/
    DynamicList<faceZone*> fz;

/*
    cell zones
    1) moving cells inside piston

    nValves*
    1) moving cells in the top of the valve
    2) moving cells in the bottom of the valve
*/

    DynamicList<cellZone*> cz;

    label nPointZones = 0;
    label nFaceZones = 0;
    label nCellZones = 0;

/*
    Adding the following faces zones:
    1:  pistonLayerFaces
    nV: pistonLayerFacesV

    Adding the following cell zones:
    1:  movingCellsPiston
    nV:  movingCellsPistonV

    Adding the following point zones:
    1: pistonPoints
    nV: valvePistonPointsV

*/

#   include "addPistonLayerFaces.H"

/*
    Adding the following face zones:

    nV: curtainCylZoneV
    nV: curtainPortZoneV
    nV: cutFaceZoneV
    nV: poppetZoneV
    nV: bottomZoneV

    Adding the following point zones:

    nV: cutPointsV
*/


#   include "addValvesFacesPointZonesPistonSliding.H"

#   include "addBowlZonesPistonSliding.H"

/*

    Adding the following point zones:

    nV: valveTopPointsV
    nV: valveBottomPointsV

    Adding the following cell zones:

    nV: movingCellsTopV
    nV: movingCellsBotV

*/

//#   include "addValvePistonCellZonesPistonSliding.H"

//#   include "addAttachDetachFacesPistonSliding.H"

    Info<< "Adding " << nPointZones << " point, "
        << nFaceZones << " face zones and "
        << nCellZones << " cell zones" << endl;

    pz.setSize(nPointZones);
    Info << "setSize pz" << endl;
    fz.setSize(nFaceZones);
    Info << "setSize fz" << endl;
    cz.setSize(nCellZones);
    Info << "setSize cz" << endl;

    addZones(pz, fz, cz);

#   include "addMeshModifiersPistonSliding.H"

    // Calculating the virtual positions of piston and valves

    setVirtualPositions();

    Info << " Write mesh " << endl;
    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    write();
    Info << " Mesh written " << endl;
}


// ************************************************************************* //
