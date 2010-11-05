/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-7 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    mixerGgiFvMesh

Description
    Simple mixer mesh using a GGI interface

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    mixerGgiFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mixerGgiFvMesh_H
#define mixerGgiFvMesh_H

#include "dynamicFvMesh.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mixerGgiFvMesh Declaration
\*---------------------------------------------------------------------------*/

class mixerGgiFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;

        // - Rotational speed in rotations per minute (rpm)
        scalar rpm_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixerGgiFvMesh(const mixerGgiFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const mixerGgiFvMesh&);


        //- Add mixer zones
        void addZonesAndModifiers();

        //- Calculate moving masks
        void calcMovingMasks() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;


public:

    //- Runtime type information
    TypeName("mixerGgiFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit mixerGgiFvMesh(const IOobject& io);


    // Destructor

        virtual ~mixerGgiFvMesh();


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Update the mesh for both mesh motion
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //