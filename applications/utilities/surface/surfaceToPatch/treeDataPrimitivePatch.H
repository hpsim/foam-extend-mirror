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
    Foam::treeDataPrimitivePatch

Description
    Encapsulation of data needed to search on PrimitivePatches

SourceFiles
    treeDataPrimitivePatch.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataPrimitivePatch_H
#define treeDataPrimitivePatch_H

#include "PrimitivePatchTemplate.H"
//#include "indexedOctree.H"
#include "treeBoundBoxList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class indexedOctree;


/*---------------------------------------------------------------------------*\
                 Class treeDataPrimitivePatchName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(treeDataPrimitivePatch);


/*---------------------------------------------------------------------------*\
                           Class treeDataPrimitivePatch Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType=point
>
class treeDataPrimitivePatch
:
    public treeDataPrimitivePatchName
{
    // Static data

        //- tolerance on linear dimensions
        static scalar tolSqr;

    // Private data

        //- Underlying geometry
        const PrimitivePatch<Face, FaceList, PointField, PointType>& patch_;

        //- Whether to precalculate and store face bounding box
        const bool cacheBb_;

        //- face bounding boxes (valid only if cacheBb_)
        treeBoundBoxList bbs_;


    // Private Member Functions

        //- Calculate face bounding box
        static treeBoundBox calcBb(const pointField&, const face&);

        //- Initialise all member data
        void update();

public:

    // Constructors

        //- Construct from patch.
        treeDataPrimitivePatch
        (
            const bool cacheBb,
            const PrimitivePatch<Face, FaceList, PointField, PointType>&
        );


    // Member Functions

        // Access

            label size() const
            {
                return patch_.size();
            }

            //- Get representative point cloud for all shapes inside
            //  (one point per shape)
            pointField points() const;


        // Search

            //- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
            //  Only makes sense for closed surfaces.
            label getVolumeType
            (
                const indexedOctree
                <
                    treeDataPrimitivePatch
                    <
                        Face,
                        FaceList,
                        PointField,
                        PointType
                    >
                >&,
                const point&
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Calculates nearest (to sample) point in shape.
            //  Returns actual point and distance (squared)
            void findNearest
            (
                const labelList& indices,
                const point& sample,

                scalar& nearestDistSqr,
                label& nearestIndex,
                point& nearestPoint
            ) const;

            //- Calculates nearest (to line) point in shape.
            //  Returns point and distance (squared)
            void findNearest
            (
                const labelList& indices,
                const linePointRef& ln,

                treeBoundBox& tightest,
                label& minIndex,
                point& linePoint,
                point& nearestPoint
            ) const
            {
                notImplemented
                (
                    "treeDataPrimitivePatch::findNearest"
                    "(const labelList&, const linePointRef&, ..)"
                );
            }

            //- Calculate intersection of shape with ray. Sets result
            //  accordingly
            bool intersects
            (
                const label index,
                const point& start,
                const point& end,
                point& result
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "treeDataPrimitivePatch.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
