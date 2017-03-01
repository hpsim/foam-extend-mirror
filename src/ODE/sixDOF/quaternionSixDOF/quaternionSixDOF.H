/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    quaternionSixDOF

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    quaternionSixDOF.C

\*---------------------------------------------------------------------------*/

#ifndef quaternionSixDOF_H
#define quaternionSixDOF_H

#include "sixDOFODE.H"
#include "finiteRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class quaternionSixDOF Declaration
\*---------------------------------------------------------------------------*/

class quaternionSixDOF
:
    public sixDOFODE
{
    // Private data

        // Body state variables

            //- Displacement relative to spring equilibrium
            dimensionedVector Xrel_;

            //- Velocity of mass centroid
            dimensionedVector U_;

            //- Average velocity of mass centroid (evaluated at midstep)
            dimensionedVector Uaverage_;

            //- Finite rotation
            finiteRotation rotation_;

            //- Rotational velocity about mass centroid
            dimensionedVector omega_;

            //- Average rotational velocity in relative coordinate system
            //  (evaluated at midstep)
            dimensionedVector omegaAverage_;

            //- Average rotational velocity in absolute coordinate system
            dimensionedVector omegaAverageAbsolute_;


        //- ODE coefficients
        scalarField coeffs_;


        // Motion constraints (given as fixed motion components)

            //- Fixed surge (x-translation)
            Switch fixedSurge_;

            //- Fixed sway (y-translation)
            Switch fixedSway_;

            //- Fixed heave (z-translation)
            Switch fixedHeave_;

            //- Fixed roll (rotation around x)
            Switch fixedRoll_;

            //- Fixed pitch (rotation around y)
            Switch fixedPitch_;

            //- Fixed yaw (rotation around z)
            Switch fixedYaw_;

            //- Constraints in referent coordinate system
            Switch referentMotionConstraints_;

            //- Rotation of referent coordinate system
            HamiltonRodriguezRot referentRotation_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        quaternionSixDOF(const quaternionSixDOF&);

        //- Disallow default bitwise assignment
        void operator=(const quaternionSixDOF&);


        // Variables in relative coordinate system (solved for)

            //- Return acceleration in relative coordinate system
            //  given current values of relative displacement and velocity
            dimensionedVector A
            (
                const dimensionedVector& xR,
                const dimensionedVector& uR,
                const HamiltonRodriguezRot& rotation
            ) const;


            //- Return rotational acceleration in relative coordinate system
            //  given current values for relative rotational velocity
            dimensionedVector OmegaDot
            (
                const HamiltonRodriguezRot& rotation,
                const dimensionedVector& omega
            ) const;

            //- Return the Euler part of moment equation
            dimensionedVector E
            (
                const dimensionedVector& omega
            ) const;

            //- Constrain rotation vector in referent or global coordinate
            //  system
            void constrainRotation(vector& vec) const;

            //- Constrain translation vector in referent or global coordinate
            //  system
            void constrainTranslation(vector& vec) const;


public:

    // Run-time type information
    TypeName("quaternionSixDOF");


    // Constructors

        //- Construct from dictionary
        quaternionSixDOF(const IOobject& io);

        //- Construct quaternionSixDOF object, changing name
        quaternionSixDOF
        (
            const word& name,
            const quaternionSixDOF& qsd
        );

        //- Return a clone, changing the name
        virtual autoPtr<sixDOFODE> clone(const word& name) const;


    // Destructor

        virtual ~quaternionSixDOF();


    // Member Functions

        // Virtual interface for 6DOF motion state

            // Variables in relative coordinate system

                //- Return displacement in translated coordinate system
                //  relative to spring equilibrium
                virtual const dimensionedVector& Xrel() const;

                //- Return rotational velocity in relative coordinate system
                virtual const dimensionedVector& omega() const;


            // Displacement and velocity in the absolute coordinate system

                //- Return position of origin in absolute coordinate system
                virtual dimensionedVector X() const;

                //- Return velocity of origin
                virtual const dimensionedVector& U() const;

                //- Return average velocity of origin (evaluated at midstep)
                virtual const dimensionedVector& Uaverage() const;


            // Non-access control for setting state variables

                //- Set ODE parameters from another ODE
                virtual void setState(const sixDOFODE&);


            // Average motion per time-step

                //- Return average rotational velocity in relative coordinate
                //  system (evaluated at midstep)
                virtual const dimensionedVector& omegaAverage() const;


            // Rotations

                //- Return rotation tensor to relative coordinate system
                virtual tensor toRelative() const;

                //- Return rotation tensor to absolute coordinate system
                virtual tensor toAbsolute() const;

                //- Return transformation tensor between new and previous
                //  rotation
                virtual const tensor& rotIncrementTensor() const;


            // ODE parameters

                //- Return number of equations
                virtual label nEqns() const
                {
                    return 13;
                }

                //- Return access to coefficients
                virtual scalarField& coeffs()
                {
                    return coeffs_;
                }

                //- Return reference to coefficients
                virtual const scalarField& coeffs() const
                {
                    return coeffs_;
                }

                //- Evaluate derivatives
                virtual void derivatives
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dydx
                ) const;

                //- Evaluate Jacobian
                virtual void jacobian
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dfdx,
                    scalarSquareMatrix& dfdy
                ) const
                {
                    notImplemented
                    (
                        "quaternionSixDOF::jacobian\n"
                        "(\n"
                        "    const scalar x,\n"
                        "    const scalarField& y,\n"
                        "    scalarField& dfdx,\n"
                        "    scalarSquareMatrix& dfdy,\n"
                        ") const"
                    );
                }

                //- Update ODE after the solution, advancing by delta
                virtual void update(const scalar delta);


        // Write controls

            //- WriteData member function required by regIOobject
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //