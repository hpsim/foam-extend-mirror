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
    Foam::curveSet

Description

SourceFiles
    curveSet.C

\*---------------------------------------------------------------------------*/

#ifndef curveSet_H
#define curveSet_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class passiveParticle;
template<class Type> class Particle;

/*---------------------------------------------------------------------------*\
                           Class curveSet Declaration
\*---------------------------------------------------------------------------*/

class curveSet
:
    public sampledSet
{
    // Private data

        //- sampling points
        List<point> sampleCoords_;


    // Private Member Functions

        //- Sample till hits boundary. Called with singleParticle at position
        //  inbetween sampleCoords_[sampleI] and sampleCoords_[sampleI+1].
        //  Returns false if end of samples reached.
        bool trackToBoundary
        (
            Particle<passiveParticle>& singleParticle,
            label& sampleI,
            DynamicList<point>& samplingPts,
            dynamicLabelList& samplingCells,
            dynamicLabelList& samplingFaces,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Samples all point in sampleCoords_
        //  samplingSegments contains segmentNo for each sample.
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
    TypeName("curve");


    // Constructors

        //- Construct from components
        curveSet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& axis,
            const List<point>& samplePoints
        );

        //- Construct from dictionary
        curveSet
        (
            const word& name,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict
        );


    // Destructor

        virtual ~curveSet();


    // Member Functions

        //- Get reference point
        virtual point getRefPoint(const List<point>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
