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

Namespace
    Foam::sixDoFRigidBodyMotionConstraints

Description
    Namespace for six DoF motion constraints


Class
    Foam::sixDoFRigidBodyMotionConstraint

Description
    Base class for defining constraints for sixDoF motions

SourceFiles
    sixDoFRigidBodyMotionConstraint.C
    newDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionConstraint_H
#define sixDoFRigidBodyMotionConstraint_H

#include "objectRegistry.H"
#include "foamTime.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sixDoFRigidBodyMotion;


/*---------------------------------------------------------------------------*\
                Class sixDoFRigidBodyMotionConstraint Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionConstraint
{

protected:

    // Protected data

        //- Constraint model specific coefficient dictionary
        dictionary sDoFRBMCCoeffs_;

        //- Solution tolerance.  Meaning depends on model, usually an
        //  absolute distance or angle.
        scalar tolerance_;

        //- Relaxation factor for solution, default to one
        scalar relaxationFactor_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionConstraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDoFRigidBodyMotionConstraint,
            dictionary,
            (const dictionary& sDoFRBMCDict),
            (sDoFRBMCDict)
        );


    // Constructors

        //- Construct from the sDoFRBMCDict dictionary and Time
        sixDoFRigidBodyMotionConstraint
        (
            const dictionary& sDoFRBMCDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionConstraint> clone() const = 0;


    // Selectors

        //- Select constructed from the sDoFRBMCDict dictionary and Time
        static autoPtr<sixDoFRigidBodyMotionConstraint> New
        (
            const dictionary& sDoFRBMCDict
        );


    // Destructor

        virtual ~sixDoFRigidBodyMotionConstraint();


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
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMCDict);

        // Access

            // Return access to sDoFRBMCCoeffs
            inline const dictionary& coeffDict() const
            {
                return sDoFRBMCCoeffs_;
            }

            //- Return access to the tolerance
            inline scalar tolerance() const
            {
                return tolerance_;
            }

            //- Return access to the relaxationFactor
            inline scalar relaxationFactor() const
            {
                return relaxationFactor_;
            }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
