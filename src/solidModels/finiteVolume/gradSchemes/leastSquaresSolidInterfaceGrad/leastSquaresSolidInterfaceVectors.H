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
    Foam::leastSquaresSolidInterfaceVectors

Description
    Least-squares gradient scheme vectors
    philipc-
    Modified to use face centre values at multi-material interfaces.

SourceFiles
    leastSquaresSolidInterfaceVectors.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresSolidInterfaceVectors_H
#define leastSquaresSolidInterfaceVectors_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                     Class leastSquaresSolidInterfaceVectors Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresSolidInterfaceVectors
:
    public MeshObject<fvMesh, leastSquaresSolidInterfaceVectors>
{
    // Private data

        //- Least-squares gradient vectors
        mutable surfaceVectorField* pVectorsPtr_;
        mutable surfaceVectorField* nVectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresSolidInterfaceVectors");


    // Constructors

        //- Construct given an fvMesh
        explicit leastSquaresSolidInterfaceVectors(const fvMesh&);


    // Destructor

        virtual ~leastSquaresSolidInterfaceVectors();


    // Member functions

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const;


        //- Update after mesh motion:
        //  Delete the least square vectors when the mesh moves
        virtual bool movePoints() const;

        //- Update after topo change:
        //  Delete the least square vectors when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
