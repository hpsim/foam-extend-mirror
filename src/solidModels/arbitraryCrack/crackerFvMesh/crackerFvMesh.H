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
    crackerFvMesh

Description
    A mesh allowing topological detachment on all internal faces of the mesh.

SourceFiles
    crackerFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef crackerFvMesh_H
#define crackerFvMesh_H

#include "topoChangerFvMesh.H"
#include "regionSplit.H"
#include "mapPolyMesh.H"
#include "faceCracker.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class crackerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class crackerFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Setup dictionary
        dictionary dict_;

        //- Last topology mapper
        autoPtr<mapPolyMesh> topoChangeMap_;

        //- Crack patch ID.  Holds faces after cracking
        polyPatchID crackPatchID_;

    // Demand-driven data

        //- Disconnected regions in the mesh
        mutable regionSplit* regionsPtr_;

        //- Number of cells in regions
        mutable labelList* nCellsInRegionPtr_;

        //- Face centres of global crack
        mutable vectorField* globalCrackFaceCentresPtr_;

        //- Face sizes of global crack
        mutable scalarField* globalCrackFaceSizesPtr_;

        //- Index of first local crack face in the global crack
        mutable label localCrackStart_;

        //- Global crack face addressing
        mutable labelList* globalCrackFaceAddressingPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        crackerFvMesh(const crackerFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const crackerFvMesh&);

        //- Add cracker zones and modifiers
        void addZonesAndModifiers();

        //- Make regions
        void makeRegions() const;

        //- Make number of cells in regions
        void makeNCellsInRegion() const;

        //- Make global crack face centres
        void makeGlobalCrackFaceCentresAndSizes() const;

        //- Make global crack face addressing
        void makeGlobalCrackFaceAddressing() const;

public:

    //- Runtime type information
    TypeName("crackerFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit crackerFvMesh(const IOobject& io);


    // Destructor

        virtual ~crackerFvMesh();


    // Member Functions

        //- Set faces to break on next update
        void setBreak
        (
            const labelList& facesToBreak,
            const boolList& faceFlip,
            const labelList& coupledFacesToBreak
        );

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Get disconnected regions
        const regionSplit& regions() const;

        //- Get number of cells in region
        label nCellsInRegion(label regI) const;

        //- Get mapper
        const mapPolyMesh& topoChangeMap() const;

        //- Get global crack face centres
        const vectorField& globalCrackFaceCentres() const;

        //- Get global crack face sizes
        const scalarField& globalCrackFaceSizes() const;

        //- Get global crack face addressing
        const labelList& globalCrackFaceAddressing() const;

        //- Get local crack face start index in global crack
        label localCrackStart() const;

        //- Get global crack size
        label globalCrackSize() const;

        //- Get global crack field
        template<class Type>
        tmp<Field<Type> > globalCrackField
        (
            const Field<Type>& localCrackField
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "crackerFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
