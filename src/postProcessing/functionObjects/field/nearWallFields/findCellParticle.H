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
    Foam::findCellParticle

Description
    Particle class that finds cells by tracking

SourceFiles
    findCellParticle.C

\*---------------------------------------------------------------------------*/

#ifndef findCellParticle_H
#define findCellParticle_H

#include "particle.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class findCellParticleCloud;

/*---------------------------------------------------------------------------*\
                     Class findCellParticle Declaration
\*---------------------------------------------------------------------------*/

class findCellParticle
:
    public particle
{
    // Private data

        //- End point to track to
        point end_;

        //- Passive data
        label data_;


public:

    friend class Cloud<findCellParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<Cloud<findCellParticle> >
    {
        labelListList& cellToData_;
        List<List<point> >& cellToEnd_;

    public:

        // Constructors

            trackingData
            (
                Cloud<findCellParticle>& cloud,
                labelListList& cellToData,
                List<List<point> >& cellToEnd
            )
            :
                particle::TrackingData<Cloud<findCellParticle> >(cloud),
                cellToData_(cellToData),
                cellToEnd_(cellToEnd)
            {}


        // Member functions

            labelListList& cellToData()
            {
                return cellToData_;
            }

            List<List<point> >& cellToEnd()
            {
                return cellToEnd_;
            }
    };


    // Constructors

        //- Construct from components
        findCellParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const point& end,
            const label data
        );

        //- Construct from Istream
        findCellParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new findCellParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<findCellParticle> operator()(Istream& is) const
            {
                return autoPtr<findCellParticle>
                (
                    new findCellParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        //- Point to track to
        const point& end() const
        {
            return end_;
        }

        //- Transported label
        label data() const
        {
            return data_;
        }


        // Tracking

            //- Track all particles to their end point
            bool move(trackingData&, const scalar);


            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch
            (
                const symmetryPlanePolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackingData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const findCellParticle&);
};


template<>
inline bool contiguous<findCellParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
