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

#include "fixedOrientation.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(fixedOrientation, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        fixedOrientation,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedOrientation::fixedOrientation
(
    const dictionary& sDoFRBMCDict
)
:
    sixDoFRigidBodyMotionConstraint(sDoFRBMCDict)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedOrientation::~fixedOrientation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionConstraints::fixedOrientation::constrain
(
    const sixDoFRigidBodyMotion& motion,
    const vector& existingConstraintForce,
    const vector& existingConstraintMoment,
    scalar deltaT,
    vector& constraintPosition,
    vector& constraintForceIncrement,
    vector& constraintMomentIncrement
) const
{
    constraintMomentIncrement = vector::zero;

    scalar maxTheta = -SMALL;

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        vector axis = vector::zero;

        axis[cmpt] = 1;

        vector refDir = vector::zero;

        refDir[(cmpt + 1) % 3] = 1;

        vector predictedDir = motion.predictedOrientation
        (
            refDir,
            existingConstraintMoment,
            deltaT
        );

        // Removing any axis component from predictedDir
        predictedDir -= (axis & predictedDir)*axis;

        scalar theta = GREAT;

        scalar magPredictedDir = mag(predictedDir);

        if (magPredictedDir > VSMALL)
        {
            predictedDir /= magPredictedDir;

            theta = acos(min(predictedDir & refDir, 1.0));

            // Recalculating axis to give correct sign to angle
            axis = (refDir ^ predictedDir);

            scalar magAxis = mag(axis);

            if (magAxis > VSMALL)
            {
                axis /= magAxis;
            }
            else
            {
                axis = vector::zero;
            }
        }

        if (theta > maxTheta)
        {
            maxTheta = theta;
        }

        constraintMomentIncrement +=
           -relaxationFactor_
           *theta*axis
           *motion.momentOfInertia()[cmpt]/sqr(deltaT);
    }

    constraintPosition = motion.centreOfMass();

    constraintForceIncrement = vector::zero;

    bool converged(mag(maxTheta) < tolerance_);

    if (motion.report())
    {
        Info<< " max angle " << maxTheta
            << " force " << constraintForceIncrement
            << " moment " << constraintMomentIncrement;

        if (converged)
        {
            Info<< " converged";
        }
        else
        {
            Info<< " not converged";
        }

        Info<< endl;
    }

    return converged;
}


bool Foam::sixDoFRigidBodyMotionConstraints::fixedOrientation::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedOrientation::write
(
    Ostream& os
) const
{
}

// ************************************************************************* //
