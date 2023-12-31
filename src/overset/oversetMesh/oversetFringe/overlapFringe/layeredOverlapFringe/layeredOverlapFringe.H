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
    layeredOverlapFringe

Description
    Layered overlap fringe algorithm. Uses neighbourhood walk nLayers deep to
    run away from user specified holes (holePatches and optional manual cell set
    containing holes) and identify acceptors.

    Holds and manages following data:

      - fringeHolesPtr_: holds fringe holes
      - acceptorsPtr_: holds acceptors
      - finalDonorAcceptorsPtr_: holds final set of donor/acceptors
      - fringeIter_: iteration counter

    No need to control the iterative process in updateIteration member function,
    acceptors are found in calcAddressing as in e.g. faceCells fringe.
    Member function updateIteration controls the iterative process:

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    layeredOverlapFringe.C

\*---------------------------------------------------------------------------*/

#ifndef layeredOverlapFringe_H
#define layeredOverlapFringe_H

#include "oversetFringe.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class layeredOverlapFringe Declaration
\*---------------------------------------------------------------------------*/

class layeredOverlapFringe
:
    public oversetFringe
{
    // Private data

        //- Fringe hole cells
        mutable labelList* fringeHolesPtr_;

        //- Acceptor cells
        mutable labelList* acceptorsPtr_;

        //- Final donor/acceptor pairs for this region (fringe)
        mutable donorAcceptorList* finalDonorAcceptorsPtr_;

        //- Number of layers
        const label nLayers_;

        //- Iteration counter
        mutable label layerIter_;


        // Optional initialization data

            //- Name of the cell set defining initial holes (empty by default).
            //  Useful when the resolution of the background mesh is much
            //  coarser than the front mesh and no hole is found
            const word holesSetName_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            layeredOverlapFringe(const layeredOverlapFringe&);

            //- Disallow default bitwise assignment
            void operator=(const layeredOverlapFringe&);


        // Addressing calculation functions

            //- Calculate hole-acceptor addressing
            void calcAddressing() const;

            //- Clear hole-acceptor addressing
            void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("layeredOverlap");


    // Static Member Functions

        //- Helper function used to evaluate coupled boundaries, excluding
        //  overset patch
        static void evaluateNonOversetBoundaries
        (
            volScalarField::GeometricBoundaryField& psib
        );


    // Constructors

        //- Construct from dictionary
        layeredOverlapFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    //- Destructor
    virtual ~layeredOverlapFringe();


    // Member Functions

        //- Update iteration.
        virtual bool updateIteration
        (
            donorAcceptorList& donorAcceptorRegionData
        ) const;

        //- Return list of hole cells
        //  Fringe hole cells are collected in addition to geometric hole
        //  cells, which fall outside of all donor regions
        virtual const labelList& fringeHoles() const;

        //- Return list of acceptor cells
        virtual const labelList& candidateAcceptors() const;

        //- Return list of final donor acceptor pairs. Note: caller may
        //  invalidate finalDonorAcceptorsPtr_ for optimisation purposes
        virtual donorAcceptorList& finalDonorAcceptors() const;

        //- Update the fringe
        virtual void update() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
