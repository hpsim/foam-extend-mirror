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
    Hamilton-Rodriguez coordinate rotation

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

Description
    Rotation defined with 4 parameters: quaternions.

\*---------------------------------------------------------------------------*/

#ifndef HamiltonRodriguezRot_H
#define HamiltonRodriguezRot_H

#include "scalarMatrices.H"
#include "coordinateRotation.H"
#include "scalarField.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class HamiltonRodriguezRot Declaration
\*---------------------------------------------------------------------------*/

class HamiltonRodriguezRot
:
    public coordinateRotation
{
    // Private data

        //- Hamilton - Rodrigez (Euler) parameters
        scalar e0_;
        scalar e1_;
        scalar e2_;
        scalar e3_;

        //- Body to quaternions transformation
        scalarRectangularMatrix Gt_;

        //- Inertial to rotated coordinate system transformation
        tensor R_;


    // Private member functions

        //- Calculate R_ - inertial to body coord. sys. rotation
        inline void calculateR()
        {
            R_.xx() = 2*(e1_*e1_ + e0_*e0_ - 0.5);
            R_.xy() = 2*(e1_*e2_ + e0_*e3_);
            R_.xz() = 2*(e1_*e3_ - e0_*e2_);

            R_.yx() = 2*(e2_*e1_ - e0_*e3_);
            R_.yy() = 2*(e2_*e2_ + e0_*e0_ - 0.5);
            R_.yz() = 2*(e2_*e3_ + e0_*e1_);

            R_.zx() = 2*(e3_*e1_ + e0_*e2_);
            R_.zy() = 2*(e3_*e2_ - e0_*e1_);
            R_.zz() = 2*(e3_*e3_ + e0_*e0_ - 0.5);
        }

        //- Calculate Gt()
        inline void calculateGt()
        {
            Gt_[0][0] = -e1();
            Gt_[0][1] = -e2();
            Gt_[0][2] = -e3();

            Gt_[1][0] =  e0();
            Gt_[1][1] = -e3();
            Gt_[1][2] =  e2();

            Gt_[2][0] =  e3();
            Gt_[2][1] =  e0();
            Gt_[2][2] = -e1();

            Gt_[3][0] = -e2();
            Gt_[3][1] =  e1();
            Gt_[3][2] =  e0();
        }


public:

    // Constructors

        //- Construct from 4 quaternions
        HamiltonRodriguezRot
        (
            const scalar e0,
            const scalar e1,
            const scalar e2,
            const scalar e3
        )
        :
            e0_(e0),
            e1_(e1),
            e2_(e2),
            e3_(e3),
            Gt_(4, 3),
            R_(tensor::zero)
        {
            calculateR();
            calculateGt();
        }

        //- Construct from rotation vector and angle
        HamiltonRodriguezRot(const vector& rVect, const scalar& rAngle)
        :
            Gt_(4, 3),
            R_(tensor::zero)
        {
            // Make rotational vector unit size
            vector r = rVect;
            scalar magR = mag(rVect);

            if (magR < SMALL)
            {
                r = vector(1, 0, 0);
            }
            else
            {
                r = rVect/magR;
            }

            e0_ = Foam::cos(0.5*rAngle);
            e1_ = Foam::sin(0.5*rAngle)*r.x();
            e2_ = Foam::sin(0.5*rAngle)*r.y();
            e3_ = Foam::sin(0.5*rAngle)*r.z();

            calculateR();
            calculateGt();
        }

        //- Construct from rotation tensor
        explicit HamiltonRodriguezRot(const tensor& R)
        :
            Gt_(4, 3),
            R_(R)
        {
            // Calculate Hamilton - Rodriguez (Euler) parameters from rotation
            // matrix

            // Note: sign of e0_ assumed positive
            e0_ = Foam::sqrt((tr(R) + 1.0)/4.0);

            // Helper variable
            const scalar oneByFourEo = 1.0/(4.0*e0_);

            e1_ = oneByFourEo*(R.zy() - R.yz());
            e2_ = oneByFourEo*(R.xz() - R.zx());
            e3_ = oneByFourEo*(R.yx() - R.xy());

            // Calculate Gt
            calculateGt();
        }


    // Destructor

        virtual ~HamiltonRodriguezRot()
        {}


    // Member functions

        // Access

            //- Returns quaternions
            inline const scalar& e0() const
            {
                return e0_;
            }

            inline const scalar& e1() const
            {
                return e1_;
            }

            inline const scalar& e2() const
            {
                return e2_;
            }

            inline const scalar& e3() const
            {
                return e3_;
            }

            //- Return initial-to-rotated transformation tensor
            virtual const tensor& R() const
            {
                return R_;
            }


            //- Return from rotated to inertial coordinate system
            inline tensor invR() const
            {
                return R().T();
            }


            inline scalar eDot
            (
                const vector& omega,
                const label quaternion
            ) const
            {
                scalar eDot = 0;

                for(label i = 0; i < 3; i++)
                {
                    eDot += Gt_[quaternion][i]*omega[i];
                }

                return 0.5*eDot;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
