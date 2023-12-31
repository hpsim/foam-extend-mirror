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

#include "thoboisSliding.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thoboisSliding::checkAndCalculate()
{

    label pistonIndex = -1;
    bool foundPiston = false;

    label linerIndex = -1;
    bool foundLiner = false;

    label cylinderHeadIndex = -1;
    bool foundCylinderHead = false;


    forAll(boundary(), i)
    {
        if (boundary()[i].name() == "piston")
        {
            pistonIndex = i;
            foundPiston = true;
        }
        else if (boundary()[i].name() == "liner")
        {
            linerIndex = i;
            foundLiner = true;
        }
        else if (boundary()[i].name() == "cylinderHead")
        {
            cylinderHeadIndex = i;
            foundCylinderHead = true;
        }
    }

    reduce(foundPiston, orOp<bool>());
    reduce(foundLiner, orOp<bool>());
    reduce(foundCylinderHead, orOp<bool>());

    if (!foundPiston)
    {
        FatalErrorIn("Foam::verticalValvesGambit::checkAndCalculate()")
            << " : cannot find piston patch"
            << abort(FatalError);
    }

    if (!foundLiner)
    {
        FatalErrorIn("Foam::verticalValvesGambit::checkAndCalculate()")
            << " : cannot find liner patch"
            << abort(FatalError);
    }

    if (!foundCylinderHead)
    {
        FatalErrorIn("Foam::verticalValvesGambit::checkAndCalculate()")
            << " : cannot find cylinderHead patch"
            << exit(FatalError);
    }

    {
        if (linerIndex != -1)
        {
            pistonPosition() =
                max(boundary()[pistonIndex].patch().localPoints()).z();
        }
        reduce(pistonPosition(), minOp<scalar>());

        if (cylinderHeadIndex != -1)
        {
            deckHeight() = min
            (
                boundary()[cylinderHeadIndex].patch().localPoints()
            ).z();

        }
        reduce(deckHeight(), minOp<scalar>());

        Info<< "deckHeight: " << deckHeight() << nl
            << "piston position: " << pistonPosition() << endl;
    }


}

void Foam::thoboisSliding::setVirtualPositions()
{
/*
    {
        virtualPistonPosition() = -GREAT;

        label pistonFaceIndex = faceZones().findZoneID("pistonLayerFaces");

        bool foundPistonFace = (pistonFaceIndex != -1);

        if(!foundPistonFace)
        {
            FatalErrorIn("Foam::verticalValvesGambit::setVirtualPistonPosition()")
                << " : cannot find the pistonLayerFaces"
                << exit(FatalError);

        }

        const labelList& pistonFaces = faceZones()[pistonFaceIndex];
        forAll(pistonFaces, i)
        {
            const face& f = faces()[pistonFaces[i]];

            // should loop over facepoints...
            forAll(f, j)
            {
                virtualPistonPosition() = max(virtualPistonPosition(), points()[f[j]].z());
            }
        }

        reduce(virtualPistonPosition(), maxOp<scalar>());

    }
*/
    forAll(valves_, valveI)
    {
        if(valves_[valveI].poppetPatchID().active())
        {

            valveTopPosition_[valveI] = -GREAT;

            label valveFaceIndex = faceZones().findZoneID("layeringFacesTopZoneV" + Foam::name(valveI + 1));

            bool foundValveFace = (valveFaceIndex != -1);

            if(!foundValveFace)
            {
                FatalErrorIn("Foam::verticalValvesGambit::setVirtualPositions()")
                    << " : cannot find the poppetValveFaces"
                    << exit(FatalError);
            }

            const labelList& valvePoppetFaces = faceZones()[valveFaceIndex];

            forAll(valvePoppetFaces, i)
            {
                const face& f = faces()[valvePoppetFaces[i]];

                // should loop over facepoints...
                forAll(f, j)
                {
                    valveTopPosition_[valveI] = max(valveTopPosition_[valveI], points()[f[j]].z());
                }
            }

            reduce(valveTopPosition_[valveI], maxOp<scalar>());


        }
        if(valves_[valveI].bottomPatchID().active())
        {
            valveBottomPosition_[valveI] = GREAT;

            label valveFaceIndex = faceZones().findZoneID("layeringFacesBottomZoneV" + Foam::name(valveI + 1));

            bool foundValveFace = (valveFaceIndex != -1);

            if(!foundValveFace)
            {
                FatalErrorIn("Foam::verticalValvesGambit::setVirtualPositions()")
                    << " : cannot find the bottomValveFaces"
                    << exit(FatalError);
            }

            const labelList& valveBottomFaces = faceZones()[valveFaceIndex];
            forAll(valveBottomFaces, i)
            {
                const face& f = faces()[valveBottomFaces[i]];

                // should loop over facepoints...
                forAll(f, j)
                {
                    valveBottomPosition_[valveI] = min(valveBottomPosition_[valveI], points()[f[j]].z());
                }
            }


        }

/*        {
            valvePistonPosition_[valveI] = -GREAT;

            label valveFaceIndex = faceZones().findZoneID("pistonLayerFacesV" + Foam::name(valveI + 1));

            const labelList& valvePistonFaces = faceZones()[valveFaceIndex];
            forAll(valvePistonFaces, i)
            {
                const face& f = faces()[valvePistonFaces[i]];

                // should loop over facepoints...
                forAll(f, j)
                {
                    valvePistonPosition_[valveI] = max(valvePistonPosition_[valveI], points()[f[j]].z());
                }
            }

            reduce(valvePistonPosition_[valveI], maxOp<scalar>());

        }
*/
    }



}
