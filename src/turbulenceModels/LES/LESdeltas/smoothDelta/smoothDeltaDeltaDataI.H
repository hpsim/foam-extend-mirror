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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Update this with w2 if applicable
inline bool smoothDelta::deltaData::update
(
    const smoothDelta::deltaData& w2,
    const scalar scale,
    const scalar tol
)
{
    if (!valid() || (delta_ < VSMALL))
    {
        // My delta not set. Take over neighbour.
        delta_ = w2.delta()/scale;

        // Something changed. Let caller know.
        return true;
    }
    else if (w2.delta() > (1 + tol)*scale*delta_)
    {
        // Neighbour is too big for me. Up my delta.
        delta_ = w2.delta()/scale;

        // Something changed. Let caller know.
        return true;
    }
    else
    {
        // Neighbour is not too big for me or change is too small
        // Nothing changed.
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline smoothDelta::deltaData::deltaData()
:
    delta_(-GREAT)
{}


// Construct from components
inline smoothDelta::deltaData::deltaData(const scalar delta)
:
    delta_(delta)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool smoothDelta::deltaData::valid() const
{
    return delta_ > -SMALL;
}


// Checks for cyclic faces
inline bool smoothDelta::deltaData::sameGeometry
(
    const polyMesh&,
    const deltaData&,
    const scalar
) const
{
    return true;
}


inline void smoothDelta::deltaData::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point&
)
{}


inline void smoothDelta::deltaData::transform
(
    const polyMesh&,
    const tensor&
)
{}


// Update absolute geometric quantities.
inline void smoothDelta::deltaData::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point&
)
{}


// Update this (cellI) with face information.
inline bool smoothDelta::deltaData::updateCell
(
    const polyMesh&,
    const label,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol
)
{
    // Take over info from face if more than deltaRatio larger.
    return update(neighbourWallInfo, maxDeltaRatio, tol);
}


// Update this (face) with cell information.
inline bool smoothDelta::deltaData::updateFace
(
    const polyMesh&,
    const label,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol
)
{
    // Take over information from cell without any scaling (scale = 1.0)
    return update(neighbourWallInfo, 1.0, tol);
}


// Update this (face) with coupled face information.
inline bool smoothDelta::deltaData::updateFace
(
    const polyMesh&,
    const label,
    const deltaData& neighbourWallInfo,
    const scalar tol
)
{
    // Take over information from coupled face without any scaling (scale = 1.0)
    return update(neighbourWallInfo, 1.0, tol);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool smoothDelta::deltaData::operator==
(
    const deltaData& rhs
) const
{
    return delta_ == rhs.delta();
}


inline bool smoothDelta::deltaData::operator!=
(
    const deltaData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
