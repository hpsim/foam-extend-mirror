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
    wedgeFaPatch

Description
    Wedge front and back plane patch.

SourceFiles
    wedgeFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFaPatch_H
#define wedgeFaPatch_H

#include "faPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wedgeFaPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgeFaPatch
:
    public faPatch
{
    // Private data

        const wedgePolyPatch* wedgePolyPatchPtr_;

        //- Axis point labels
        mutable labelList axisPoints_;

        //- Have axis points already been checked?
        mutable bool axisPointsChecked_;

        //- Find axis points
        void findAxisPoints() const;

public:

    //- Runtime type information
    TypeName("wedge");


    // Constructors

        //- Construct from dictionary
        wedgeFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );


    // Destructor

        virtual ~wedgeFaPatch()
        {}


    // Member functions

        // Access

            //- Return axis of the wedge
            const vector& axis() const
            {
                return wedgePolyPatchPtr_->axis();
            }

            //- Return plane normal between the wedge boundaries
            const vector& centreNormal() const
            {
                return wedgePolyPatchPtr_->centreNormal();
            }

            //- Return face transformation tensor
            const tensor& edgeT() const
            {
                return wedgePolyPatchPtr_->faceT();
            }

            //- Return neighbour-cell transformation tensor
            const tensor& faceT() const
            {
                return wedgePolyPatchPtr_->cellT();
            }

            //- Return first axis point label
            label axisPoint() const
            {
                if (!axisPointsChecked_)
                {
                    findAxisPoints();
                }

                if (axisPoints_.size() > 0)
                {
                    return axisPoints_[0];
                }
                else
                {
                    return -1;
                }
            }

            //- Return all axis point labels
            const labelList& axisPoints() const
            {
                if (!axisPointsChecked_)
                {
                    findAxisPoints();
                }

                return axisPoints_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
