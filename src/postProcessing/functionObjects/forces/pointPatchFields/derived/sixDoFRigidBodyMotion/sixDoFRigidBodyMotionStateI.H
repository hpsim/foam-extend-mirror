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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const Foam::point& Foam::sixDoFRigidBodyMotionState::centreOfMass() const
{
    return centreOfMass_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotionState::Q() const
{
    return Q_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionState::v() const
{
    return v_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionState::a() const
{
    return a_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionState::pi() const
{
    return pi_;
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotionState::tau() const
{
    return tau_;
}


inline Foam::point& Foam::sixDoFRigidBodyMotionState::centreOfMass()
{
    return centreOfMass_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotionState::Q()
{
    return Q_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionState::v()
{
    return v_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionState::a()
{
    return a_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionState::pi()
{
    return pi_;
}


inline Foam::vector& Foam::sixDoFRigidBodyMotionState::tau()
{
    return tau_;
}


// ************************************************************************* //
