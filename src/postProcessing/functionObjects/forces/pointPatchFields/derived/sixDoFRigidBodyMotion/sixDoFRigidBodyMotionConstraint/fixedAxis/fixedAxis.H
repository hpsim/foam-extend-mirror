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
    Foam::sixDoFRigidBodyMotionConstraints::fixedAxis

Description
    sixDoFRigidBodyMotionConstraint.  Axis of body fixed global
    space.

SourceFiles
    fixedAxis.C

\*---------------------------------------------------------------------------*/

#ifndef fixedAxis_H
#define fixedAxis_H

#include "sixDoFRigidBodyMotionConstraint.H"
#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionConstraints
{

/*---------------------------------------------------------------------------*\
                            Class fixedAxis Declaration
\*---------------------------------------------------------------------------*/

class fixedAxis
:
    public sixDoFRigidBodyMotionConstraint
{

    // Private data

        //- Reference axis in global space
        vector fixedAxis_;


public:

    //- Runtime type information
    TypeName("fixedAxis");


    // Constructors

        //- Construct from components
        fixedAxis
        (
            const dictionary& sDoFRBMCDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionConstraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionConstraint>
            (
                new fixedAxis(*this)
            );
        }


    // Destructor

        virtual ~fixedAxis();


    // Member Functions

        //- Calculate the constraint position, force and moment.
        //  Global reference frame vectors.  Returns boolean stating
        //  whether the constraint been converged to tolerance.
        virtual bool constrain
        (
            const sixDoFRigidBodyMotion& motion,
            const vector& existingConstraintForce,
            const vector& existingConstraintMoment,
            scalar deltaT,
            vector& constraintPosition,
            vector& constraintForceIncrement,
            vector& constraintMomentIncrement
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
