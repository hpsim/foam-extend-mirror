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
    dieselEngineValve

Description

SourceFiles
    dieselEngineValve.C

\*---------------------------------------------------------------------------*/

#ifndef dieselEngineValve_H
#define dieselEngineValve_H

#include "word.H"
#include "coordinateSystem.H"
#include "polyPatchID.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                          Class dieselEngineValve Declaration
\*---------------------------------------------------------------------------*/

class dieselEngineValve
{
    // Private data

        //- Name of valve
        word name_;

        //- Reference to engine mesh
        const polyMesh& mesh_;

        //- Reference to engine database
        const engineTime& engineDB_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Patch and zone names

            //- Valve bottom patch
            polyPatchID bottomPatch_;

            //- Valve poppet patch
            polyPatchID poppetPatch_;

            //- Valve poppet patch
            polyPatchID sidePatch_;

            //- Valve stem patch
            polyPatchID stemPatch_;

            //- Valve curtain manifold patch
            polyPatchID downInPortPatch_;

            //- Valve curtain cylinder patch
            polyPatchID downInCylinderPatch_;

            //- Valve curtain manifold patch
            polyPatchID upInPortPatch_;

            //- Valve curtain cylinder patch
            polyPatchID upInCylinderPatch_;

            //- Valve detach in cylinder patch
            polyPatchID detachInCylinderPatch_;

            //- Valve detach in port patch
            polyPatchID detachInPortPatch_;

            //- Faces to detach
            labelList detachFaces_;

            //- detach tolerance
            scalar detachTol_;

        // Valve lift data

            //- Valve lift profile
            graph liftProfile_;

            //- Lift curve start angle
            scalar liftProfileStart_;

            //- Lift curve end angle
            scalar liftProfileEnd_;

            //- Minimum valve lift.  On this lift the valve is considered closed
            const scalar minLift_;


            //- Valve diameter
            const scalar diameter_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        dieselEngineValve(const dieselEngineValve&);

        //- Disallow default bitwise assignment
        void operator=(const dieselEngineValve&);


        //- Adjust crank angle to drop within the limits of the lift profile
        scalar adjustCrankAngle(const scalar theta) const;

public:

    // Static data members


    // Constructors

        //- Construct from components
        dieselEngineValve
        (
            const word& name,
            const polyMesh& mesh,
            const autoPtr<coordinateSystem>& valveCS,
            const word& bottomPatchName,
            const word& poppetPatchName,
            const word& sidePatchName,
            const word& stemPatchName,
            const word& downInPortPatchName,
            const word& downInCylinderPatchName,
            const word& upInPortPatchName,
            const word& upInCylinderPatchName,
            const word& detachInCylinderPatchName,
            const word& detachInPortPatchName,
            const labelList& detachFaces,
            const scalar& detachTol,
            const graph& liftProfile,
            const scalar minLift,
            const scalar diameter

        );

        //- Construct from dictionary
        dieselEngineValve
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return lift profile
        const graph& liftProfile() const
        {
            return liftProfile_;
        }

        //- Return valve diameter
        scalar diameter() const
        {
            return diameter_;
        }


        // Valve patches

            //- Return ID of bottom patch
            const polyPatchID& bottomPatchID() const
            {
                return bottomPatch_;
            }

            //- Return ID of poppet patch
            const polyPatchID& poppetPatchID() const
            {
                return poppetPatch_;
            }

            //- Return ID of side patch
            const polyPatchID& sidePatchID() const
            {
                return sidePatch_;
            }

            //- Return ID of stem patch
            const polyPatchID& stemPatchID() const
            {
                return stemPatch_;
            }

            //- Return ID of down surface in cylinder patch
            const polyPatchID& downInCylinderPatchID() const
            {
                return downInCylinderPatch_;
            }

            //- Return ID of down surface  in port patch
            const polyPatchID& downInPortPatchID() const
            {
                return downInPortPatch_;
            }

            //- Return ID of up surface  in cylinder patch
            const polyPatchID& upInCylinderPatchID() const
            {
                return upInCylinderPatch_;
            }

            //- Return ID of up surface  in port patch
            const polyPatchID& upInPortPatchID() const
            {
                return upInPortPatch_;
            }

            //- Return ID of detach in cylinder patch
            const polyPatchID& detachInCylinderPatchID() const
            {
                return detachInCylinderPatch_;
            }

            //- Return ID of detach in port patch
            const polyPatchID& detachInPortPatchID() const
            {
                return detachInPortPatch_;
            }

            //- Return face labels of detach curtain
            const labelList& detachFaces() const
            {
                return detachFaces_;
            }

            //- Return face labels of detach curtain
            const scalar& detachTol() const
            {
                return detachTol_;
            }

        // Valve position and velocity

            //- Return valve lift given crank angle in degrees
            scalar lift(const scalar theta) const;

            //- Is the valve open?
            bool isOpen() const;

            //- Return current lift
            scalar curLift() const;

            //- Return valve velocity for current time-step
            scalar curVelocity() const;

            //- Return list of active patch labels for the valve head
            //  (stem is excluded)
            labelList movingPatchIDs() const;


        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
