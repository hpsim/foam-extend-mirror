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

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label referredCell::sourceProc() const
{
    return sourceProc_;
}


inline label referredCell::sourceCell() const
{
    return sourceCell_;
}


inline const vector& referredCell::offset() const
{
    return offset_;
}


inline const tensor& referredCell::rotation() const
{
    return rotation_;
}


inline const vectorList& referredCell::vertexPositions() const
{
    return vertexPositions_;
}


inline const edgeList& referredCell::edges() const
{
    return edges_;
}


inline const labelListList& referredCell::faces() const
{
    return faces_;
}


inline const vectorList& referredCell::faceCentres() const
{
    return faceCentres_;
}


inline const vectorList& referredCell::faceAreas() const
{
    return faceAreas_;
}


inline labelList& referredCell::realCells()
{
    return realCellsForInteraction_;
}


inline const labelList& referredCell::realCellsForInteraction() const
{
    return realCellsForInteraction_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool operator==
(
    const referredCell& a,
    const referredCell& b
)
{
    return const_cast<referredCell&>(a).duplicate
    (
        const_cast<const referredCell&>(b)
    );
}


inline bool operator!=
(
    const referredCell& a,
    const referredCell& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
