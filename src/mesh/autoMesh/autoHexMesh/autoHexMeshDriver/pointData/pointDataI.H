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

#include "polyMesh.H"
#include "transform.H"
#include "wallPoint.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
inline bool Foam::pointData::update
(
    const point& pt,
    const pointData& w2,
    const scalar tol
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    if (!valid())
    {
        distSqr_ = dist2;
        origin_ = w2.origin();
        s_ = w2.s();
        v_ = w2.v();

        return true;
    }


//    if (v_ != w2.v())
//    {
//        return false;
//    }


    scalar diff = distSqr_ - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr_ > SMALL) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ = dist2;
        origin_ = w2.origin();
        s_ = w2.s();
        v_ = w2.v();

        return true;
    }
}


// Update this with w2 (information on same point)
inline bool Foam::pointData::update
(
    const pointData& w2,
    const scalar tol
)
{
    if (!valid())
    {
        // current not yet set so use any value
        distSqr_ = w2.distSqr();
        origin_ = w2.origin();
        s_ = w2.s();
        v_ = w2.v();
        return true;
    }


//    if (v_ != w2.v())
//    {
//        return false;
//    }


    scalar diff = distSqr_ - w2.distSqr();

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr_ > SMALL) && (diff/distSqr_ < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // update with new values
        distSqr_ =  w2.distSqr();
        origin_ = w2.origin();
        s_ = w2.s();
        v_ = w2.v();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::pointData::pointData()
:
    origin_(wallPoint::greatPoint),
    distSqr_(GREAT),
    s_(GREAT),
    v_(wallPoint::greatPoint)
{}


// Construct from origin, distance
inline Foam::pointData::pointData
(
    const point& origin,
    const scalar distSqr,
    const scalar s,
    const vector& v
)
:
    origin_(origin),
    distSqr_(distSqr),
    s_(s),
    v_(v)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point& Foam::pointData::origin() const
{
    return origin_;
}


inline Foam::scalar Foam::pointData::distSqr() const
{
    return distSqr_;
}


inline Foam::scalar Foam::pointData::s() const
{
    return s_;
}


inline const Foam::vector& Foam::pointData::v() const
{
    return v_;
}


inline bool Foam::pointData::valid() const
{
    return origin_ != wallPoint::greatPoint;
}


// Checks for cyclic points
inline bool Foam::pointData::sameGeometry
(
    const pointData& w2,
    const scalar tol
) const
{
    scalar diff = Foam::mag(distSqr() - w2.distSqr());

    if (diff < SMALL)
    {
        return true;
    }
    else
    {
        if ((distSqr() > SMALL) && ((diff/distSqr()) < tol))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


inline void Foam::pointData::leaveDomain
(
    const polyPatch& patch,
    const label patchPointI,
    const point& coord
)
{
    origin_ -= coord;
}


inline void Foam::pointData::transform(const tensor& rotTensor)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


// Update absolute geometric quantities. Note that distance (distSqr_)
// is not affected by leaving/entering domain.
inline void Foam::pointData::enterDomain
(
    const polyPatch& patch,
    const label patchPointI,
    const point& coord
)
{
    // back to absolute form
    origin_ += coord;
}


// Update this with information from connected edge
inline bool Foam::pointData::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const label edgeI,
    const pointData& edgeInfo,
    const scalar tol
)
{
    return
        update
        (
            mesh.points()[pointI],
            edgeInfo,
            tol
        );
    }


// Update this with new information on same point
inline bool Foam::pointData::updatePoint
(
    const polyMesh& mesh,
    const label pointI,
    const pointData& newPointInfo,
    const scalar tol
)
{
    return
        update
        (
            mesh.points()[pointI],
            newPointInfo,
            tol
        );
}


// Update this with new information on same point. No extra information.
inline bool Foam::pointData::updatePoint
(
    const pointData& newPointInfo,
    const scalar tol
)
{
    return update(newPointInfo, tol);
}


// Update this with information from connected point
inline bool Foam::pointData::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointI,
    const pointData& pointInfo,
    const scalar tol
)
{
    const pointField& points = mesh.points();

    const edge& e = mesh.edges()[edgeI];

    const point edgeMid(0.5*(points[e[0]] + points[e[1]]));

    return
        update
        (
            edgeMid,
            pointInfo,
            tol
        );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointData::operator==(const pointData& rhs) const
{
    return origin() == rhs.origin();
}


inline bool Foam::pointData::operator!=(const pointData& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
