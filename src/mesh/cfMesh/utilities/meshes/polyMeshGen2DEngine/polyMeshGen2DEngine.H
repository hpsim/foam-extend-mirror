/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
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
    polyMeshGen2DEngine

Description
    A simple engine which provides topological information of a 2D mesh
    and allows for maintaining consistency

SourceFiles
    polyMeshGen2DEngine.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGen2DEngine_H
#define polyMeshGen2DEngine_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen.H"
#include "boolList.H"
#include "labelList.H"
#include "boundBox.H"
#include "demandDrivenData.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGen2DEngine
{
    // Private data
        //- const reference to the mesh
        const polyMeshGen& mesh_;

        //- bounding box
        boundBox bb_;

        //- which faces are not in the x-y plane
        mutable boolList* activeFacePtr_;

        //- labels of active faces
        mutable labelList* activeFaceLabelsPtr_;

        //- which points are in the x-y plane and at smallest z
        mutable boolList* zMinPointPtr_;

        //- labels of points at minimum z
        mutable labelList* zMinPointLabelsPtr_;

        //- labels of offset points of zMin points
        mutable labelList* zMinToZMaxPtr_;

        //- which points are in the x-y plane and at largest z
        mutable boolList* zMaxPointPtr_;

        //- labels of points at maximum z
        mutable labelList* zMaxPointLabelsPtr_;

        //- labels of zMin points of zMax points
        mutable labelList* zMaxToZMinPtr_;

    // Private member functions
        //- find active faces
        void findActiveFaces() const;

        //- find active face labels
        void findActiveFaceLabels() const;

        //- find points at minimum z
        void findZMinPoints() const;

        //- find labels of points at minimum z
        void findZMinPointLabels() const;

        //- find offset point to each zMin point
        void findZMinOffsetPoints() const;

        //- find points at maximum z
        void findZMaxPoints() const;

        //- find labels of points at minimum z
        void findZMaxPointLabels() const;

        //- find offset point to each zMin point
        void findZMaxOffsetPoints() const;

        //- disallow bitwise assignment
        void operator=(const polyMeshGen2DEngine&);

        //- disallow copy construct
        polyMeshGen2DEngine(const polyMeshGen2DEngine&);

public:

    // Constructors
        //- Construct from the const reference to the mesh
        polyMeshGen2DEngine(const polyMeshGen& mesh);

    // Destructor
        ~polyMeshGen2DEngine();

    // Member functions
        //- const access to active faces
        inline const boolList& activeFace() const;

        //- labels of active faces
        inline const labelList& activeFaceLabels() const;

        //- const access to an array of booleans which is true for points
        //- in the x-y plane with the smallest z coordinate
        inline const boolList& zMinPoints() const;

        //- labels of points in the x-y with the smallest z coordinate
        inline const labelList& zMinPointLabels() const;

        //- offset points of points at minimum z
        inline const labelList& zMinToZMax() const;

        //- const access to an array of booleans which is true for points
        //- in the x-y plane with the largest z coordinate
        inline const boolList& zMaxPoints() const;

        //- labels of points in the x-y with the largest z coordinate
        inline const labelList& zMaxPointLabels() const;

        //- offset points of points at maximum z
        inline const labelList& zMaxToZMin() const;

        //- correct x and y coordinates of offset points
        //- and unify the z coordinate
        void correctPoints();

        //- delete all dynamically allocated data
        void clearOut();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen2DEngineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
