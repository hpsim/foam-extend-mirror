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
    Foam::sixDoFRigidBodyMotion

Description
    Six degree of freedom motion for a rigid body.  Angular momentum
    stored in body fixed reference frame.  Reference orientation of
    the body (where Q = I) must align with the cartesian axes such
    that the Inertia tensor is in principle component form.

    Symplectic motion as per:

    title = {Symplectic splitting methods for rigid body molecular dynamics},
    publisher = {AIP},
    year = {1997},
    journal = {The Journal of Chemical Physics},
    volume = {107},
    number = {15},
    pages = {5840-5851},
    url = {http://link.aip.org/link/?JCP/107/5840/1},
    doi = {10.1063/1.474310}

    Can add restraints (i.e. a spring) and constraints (i.e. motion
    may only be on a plane).

SourceFiles
    sixDoFRigidBodyMotionI.H
    sixDoFRigidBodyMotion.C
    sixDoFRigidBodyMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotion_H
#define sixDoFRigidBodyMotion_H

#include "sixDoFRigidBodyMotionState.H"
#include "pointField.H"
#include "sixDoFRigidBodyMotionRestraint.H"
#include "sixDoFRigidBodyMotionConstraint.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class sixDoFRigidBodyMotion;
Istream& operator>>(Istream&, sixDoFRigidBodyMotion&);
Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotion&);


/*---------------------------------------------------------------------------*\
                      Class sixDoFRigidBodyMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotion
{
    // Private data

        //- Motion state data object
        sixDoFRigidBodyMotionState motionState_;

        //- Restraints on the motion
        PtrList<sixDoFRigidBodyMotionRestraint> restraints_;

        //- Names of the restraints
        wordList restraintNames_;

        //- Constaints on the motion
        PtrList<sixDoFRigidBodyMotionConstraint> constraints_;

        //- Names of the constraints
        wordList constraintNames_;

        //- Maximum number of iterations allowed to attempt to obey
        //  constraints
        label maxConstraintIterations_;

        //- Centre of mass of initial state
        point initialCentreOfMass_;

        //- Orientation of initial state
        tensor initialQ_;

        //- Moment of inertia of the body in reference configuration
        //  (Q = I)
        diagTensor momentOfInertia_;

        //- Mass of the body
        scalar mass_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;


    // Private Member Functions

        //- Calculate the rotation tensor around the body reference
        //  frame x-axis by the given angle
        inline tensor rotationTensorX(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame y-axis by the given angle
        inline tensor rotationTensorY(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame z-axis by the given angle
        inline tensor rotationTensorZ(scalar deltaT) const;

        //- Apply rotation tensors to Q for the given torque (pi) and deltaT
        inline void rotate(tensor& Q, vector& pi, scalar deltaT) const;

        //- Apply the restraints to the object
        void applyRestraints();

        //- Apply the constraints to the object
        void applyConstraints(scalar deltaT);

        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return access to the motion state
            inline const sixDoFRigidBodyMotionState& motionState() const;

            //- Return access to the restraints
            inline const PtrList<sixDoFRigidBodyMotionRestraint>&
                restraints() const;

            //- Return access to the restraintNames
            inline const wordList& restraintNames() const;

            //- Return access to the constraints
            inline const PtrList<sixDoFRigidBodyMotionConstraint>&
                constraints() const;

            //- Return access to the constraintNames
            inline const wordList& constraintNames() const;

            //- Return access to the maximum allowed number of
            //  constraint iterations
            inline label maxConstraintIterations() const;

            //- Return access to the initial centre of mass
            inline const point& initialCentreOfMass() const;

            //- Return access to the initial orientation
            inline const tensor& initialQ() const;

            //- Return access to the orientation
            inline const tensor& Q() const;

            //- Return access to velocity
            inline const vector& v() const;

            //- Return access to acceleration
            inline const vector& a() const;

            //- Return access to angular momentum
            inline const vector& pi() const;

            //- Return access to torque
            inline const vector& tau() const;


        // Edit

            //- Return access to the centre of mass
            inline point& initialCentreOfMass();

            //- Return access to the centre of mass
            inline tensor& initialQ();

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to vector
            inline vector& v();

            //- Return non-const access to acceleration
            inline vector& a();

            //- Return non-const access to angular momentum
            inline vector& pi();

            //- Return non-const access to torque
            inline vector& tau();


public:

    // Constructors

        //- Construct null
        sixDoFRigidBodyMotion();

        //- Construct from components
        sixDoFRigidBodyMotion
        (
            const point& centreOfMass,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            scalar mass,
            const point& initialCentreOfMass,
            const tensor& initialQ,
            const diagTensor& momentOfInertia,
            bool report = false
        );

        //- Construct from dictionary
        sixDoFRigidBodyMotion(const dictionary& dict);

        //- Construct as copy
        sixDoFRigidBodyMotion(const sixDoFRigidBodyMotion&);


    //- Destructor
    ~sixDoFRigidBodyMotion();


    // Member Functions

        //- Add restraints to the motion, public to allow external
        //  addition of restraints after construction
        void addRestraints(const dictionary& dict);

        //- Add restraints to the motion, public to allow external
        //  addition of restraints after construction
        void addConstraints(const dictionary& dict);

        //- First leapfrog velocity adjust and motion part, required
        // before force calculation
        void updatePosition
        (
            scalar deltaT
        );

        //- Second leapfrog velocity adjust part, required after motion and
        // force calculation
        void updateForce
        (
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT
        );

        //- Global forces supplied at locations, calculating net force
        //  and moment
        void updateForce
        (
            const pointField& positions,
            const vectorField& forces,
            scalar deltaT
        );

        //- Transform the given initial state pointField by the current
        //  motion state
        inline tmp<pointField> currentPosition
        (
            const pointField& pInitial
        ) const;

        //- Transform the given initial state point by the current motion
        //  state
        inline point currentPosition(const point& pInitial) const;

        //- Transform the given initial state direction by the current
        //  motion state
        inline vector currentOrientation(const vector& vInitial) const;

        //- Access the orientation tensor, Q.
        //  globalVector = Q & bodyLocalVector
        //  bodyLocalVector = Q.T() & globalVector
        inline const tensor& orientation() const;

        //- Predict the position of the supplied initial state point
        //  after deltaT given the current motion state and the
        //  additional supplied force and moment
        point predictedPosition
        (
            const point& pInitial,
            const vector& deltaForce,
            const vector& deltaMoment,
            scalar deltaT
        ) const;

        //- Predict the orientation of the supplied initial state
        //  vector after deltaT given the current motion state and the
        //  additional supplied moment
        vector predictedOrientation
        (
            const vector& vInitial,
            const vector& deltaMoment,
            scalar deltaT
        ) const;

        //- Return the angular velocity in the global frame
        inline vector omega() const;

        //- Return the velocity of a position given by the current
        //  motion state
        inline point currentVelocity(const point& pt) const;

        //- Report the status of the motion
        void status() const;


        // Access

            //- Return const access to the centre of mass
            inline const point& centreOfMass() const;

            //- Return access to the inertia tensor
            inline const diagTensor& momentOfInertia() const;

            //- Return const access to the mass
            inline scalar mass() const;

            //- Return the report Switch
            inline bool report() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& centreOfMass();

            //- Return non-const access to the inertia tensor
            inline diagTensor& momentOfInertia();

            //- Return non-const access to the mass
            inline scalar& mass();


        //- Write
        void write(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, sixDoFRigidBodyMotion&);
        friend Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotion&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
