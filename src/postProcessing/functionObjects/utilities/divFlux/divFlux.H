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
    divFlux

Description
    FunctionObject calculates the divergence of flux field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    divFlux.C

\*---------------------------------------------------------------------------*/

#ifndef divFlux_H
#define divFlux_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class divFlux Declaration
\*---------------------------------------------------------------------------*/

class divFlux
:
    public functionObject
{
    // Private data

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Flux field name
        word phiName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        divFlux
        (
            const divFlux&
        );

        //- Disallow default bitwise assignment
        void operator=(const divFlux&);


public:

    //- Runtime type information
    TypeName("divFlux");


    // Constructors

        //- Construct from components
        divFlux
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for mesh point-motion
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
