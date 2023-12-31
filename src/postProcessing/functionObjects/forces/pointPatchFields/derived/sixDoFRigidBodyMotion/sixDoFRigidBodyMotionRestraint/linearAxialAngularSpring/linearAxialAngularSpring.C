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

#include "linearAxialAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(linearAxialAngularSpring, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        linearAxialAngularSpring,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::
linearAxialAngularSpring
(
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(sDoFRBMRDict),
    refQ_(),
    axis_(),
    stiffness_(),
    damping_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::
~linearAxialAngularSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    vector refDir = rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 1, 0);

    vector oldDir = refQ_ & refDir;

    vector newDir = motion.orientation() & refDir;

    if (mag(oldDir & axis_) > 0.95 || mag(newDir & axis_) > 0.95)
    {
        // Directions getting close to the axis, change reference

        refDir = rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 0, 1);

        // Must update existing variables - not create new ones
        // HR 18/Jul/2013

        oldDir = refQ_ & refDir;

        newDir = motion.orientation() & refDir;
    }

    // Removing any axis component from oldDir and newDir and normalising
    oldDir -= (axis_ & oldDir)*axis_;
    oldDir /= (mag(oldDir) + VSMALL);

    newDir -= (axis_ & newDir)*axis_;
    newDir /= (mag(newDir) + VSMALL);

    scalar theta = mag(acos(min(oldDir & newDir, 1.0)));

    // Temporary axis with sign information.
    vector a = (oldDir ^ newDir);

    // Remove any component that is not along axis that may creep in
    a = (a & axis_)*axis_;

    scalar magA = mag(a);

    if (magA > VSMALL)
    {
        a /= magA;
    }
    else
    {
        a = vector::zero;
    }

    // Damping of along axis angular velocity only
    restraintMoment = -stiffness_*theta*a - damping_*(motion.omega() & a)*a;

    restraintForce = vector::zero;

    // Not needed to be altered as restraintForce is zero, but set to
    // centreOfMass to be sure of no spurious moment
    restraintPosition = motion.centreOfMass();

    if (motion.report())
    {
        Info<< " angle " << theta
            << " force " << restraintForce
            << " moment " << restraintMoment
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    refQ_ = sDoFRBMRCoeffs_.lookupOrDefault<tensor>("referenceOrientation", I);

    if (mag(mag(refQ_) - sqrt(3.0)) > 1e-9)
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::"
            "read"
            "("
                "const dictionary& sDoFRBMRDict"
            ")"
        )
            << "referenceOrientation " << refQ_ << " is not a rotation tensor. "
            << "mag(referenceOrientation) - sqrt(3) = "
            << mag(refQ_) - sqrt(3.0) << nl
            << exit(FatalError);
    }

    axis_ = sDoFRBMRCoeffs_.lookup("axis");

    scalar magAxis(mag(axis_));

    if (magAxis > VSMALL)
    {
        axis_ /= magAxis;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::"
            "read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "axis has zero length"
            << abort(FatalError);
    }

    sDoFRBMRCoeffs_.lookup("stiffness") >> stiffness_;

    sDoFRBMRCoeffs_.lookup("damping") >> damping_;

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring::write
(
    Ostream& os
) const
{
    os.writeKeyword("referenceOrientation")
        << refQ_ << token::END_STATEMENT << nl;

    os.writeKeyword("axis")
        << axis_ << token::END_STATEMENT << nl;

    os.writeKeyword("stiffness")
        << stiffness_ << token::END_STATEMENT << nl;

    os.writeKeyword("damping")
        << damping_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
