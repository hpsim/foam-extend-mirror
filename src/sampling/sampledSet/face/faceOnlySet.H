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
    Foam::faceOnlySet

Description

SourceFiles
    faceOnlySet.C

\*---------------------------------------------------------------------------*/

#ifndef faceOnlySet_H
#define faceOnlySet_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class passiveParticle;
template<class Type> class Particle;

/*---------------------------------------------------------------------------*\
                           Class faceOnlySet Declaration
\*---------------------------------------------------------------------------*/

class faceOnlySet
:
    public sampledSet
{
    // Private data

        //- Starting point
        point start_;

        //- End point
        point end_;


    // Static data

        //- Maximum number of steps for intersection search. Introduced to avoid
        //  inifinite loops (IG 4/Dec/2018)
        static const debug::optimisationSwitch maxNSteps_;


    // Private Member Functions

        //- Samples from startTrackPt/CellI. Updates particle/samplePt/sampleI
        //  and puts
        //  samples in the DynamicLists. Returns false if end of all samples
        //  reached
        bool trackToBoundary
        (
            Particle<passiveParticle>& singleParticle,
            DynamicList<point>& samplingPts,
            dynamicLabelList& samplingCells,
            dynamicLabelList& samplingFaces,
            DynamicList<scalar>& samplingCurve
        ) const;

        //- Samples from start_ to end_. samplingSegments contains segmentNo
        //  for each sample.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            dynamicLabelList& samplingCells,
            dynamicLabelList& samplingFaces,
            dynamicLabelList& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();

public:

    //- Runtime type information
    TypeName("face");

    // Constructors

        //- Construct from components
        faceOnlySet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& axis,
            const point& start,
            const point& end
        );

        //- Construct from dictionary
        faceOnlySet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict
        );


    // Destructor

        virtual ~faceOnlySet();


    // Member Functions

        const point& start() const
        {
            return start_;
        }

        const point& end() const
        {
            return end_;
        }

        //- Get reference point
        virtual point getRefPoint(const List<point>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
