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
    Foam::mergePolyMesh

Description
    Add a given mesh to the original mesh to create a single new mesh

SourceFiles
    mergePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mergePolyMesh_H
#define mergePolyMesh_H

#include "polyMesh.H"
#include "directTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class mergePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mergePolyMesh
:
    public polyMesh
{
    // Private data

        //- Topological change to accumulated all mesh changes
        directTopoChange meshMod_;

        //- Patch types
        DynamicList<word> patchTypes_;

        //- Patch names
        DynamicList<word> patchNames_;

        //- Point zone names
        DynamicList<word> pointZoneNames_;

        //- Face zone names
        DynamicList<word> faceZoneNames_;

        //- Cell zone names
        DynamicList<word> cellZoneNames_;

        //- Point zones
        dynamicLabelList  pointZones_;

        //- Face zones
        dynamicLabelList  faceZones_;

        //- Face zone Flips
        DynamicList<bool>  faceZoneFlips_;

        //- Cell zones
        dynamicLabelList cellZones_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mergePolyMesh(const mergePolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const mergePolyMesh&);


        //- Return patch index given a name and type
        label patchIndex(const polyPatch&);

        //- Return zone index given a list of active zones and a name
        label zoneIndex(DynamicList<word>&, const word&);


public:

    //- Runtime type information
    TypeName("mergePolyMesh");


    // Constructors

        //- Construct from IOobject
        mergePolyMesh(const IOobject& io);


    //- Destructor
    virtual ~mergePolyMesh()
    {}


    // Member Functions

        //- Add a mesh
        void addMesh(const polyMesh& m);

        //- Merge meshes
        void merge();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
