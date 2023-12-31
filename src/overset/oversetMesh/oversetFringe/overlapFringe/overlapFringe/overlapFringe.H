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
    overlapFringe

Description
    Overlap fringe algorithm

    Holds and manages following data:

      - fringeHolesPtr_: holds fringe holes which are basically discarded
        acceptors from previous iterations
      - acceptorsPtr_: holds current set of acceptors to consider
      - finalDonorAcceptorsPtr_: final set of donor/acceptors which will be used
        for the overset interpolation
      - donorSuitability_: run-time selectable Donor Suitability Function
      - initPatchNames_: optional list of patches that we use to start the
        acceptor search from
      - minGlobalFraction_: minimum allowable global fraction of suitable pairs
      - cumulativeDonorAcceptorsPtr_: list of suitable donor/acceptor pairs
      - fringeIter_: iteration counter

    Member function updateIteration controls the iterative process:
      - Checks whether given donor/acceptors are suitable and filters them,
      - If the user-specified criterion controlled by minGlobalFraction_ is not
        reached (i.e. suitable overlap has not been found), then we loop through
        unsuitable acceptors and mark their neighbours (excluding existing
        acceptors and holes) as candidate acceptors for the next iteration,
      - Unsuitable acceptors become fringe holes.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    overlapFringe.C

\*---------------------------------------------------------------------------*/

#ifndef overlapFringe_H
#define overlapFringe_H

#include "oversetFringe.H"
#include "donorSuitability.H"
#include "optimisationSwitch.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class overlapFringe Declaration
\*---------------------------------------------------------------------------*/

class overlapFringe
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


        // Optional initialization data

            //- Name of the cell set defining initial holes (empty by default).
            //  Useful when the resolution of the background mesh is much
            //  coarser than the front mesh and no hole is found
            const word holesZoneName_;

            //- Optional list of patches to start the iterative fringe assembly
            //  process (empty list by default). Useful when we actually have
            //  a hole as a patch in the background mesh
            const wordList initPatchNames_;


        // Data for iteration control

            //- Run-time selectable donor/acceptor suitability handler
            autoPtr<donorSuitability::donorSuitability> donorSuitability_;

            //- Minimum global fraction of suitable donors. Iteration process
            //  stops if the calculated fraction is greater than the minimum
            //  allowable fraction
            const scalar minGlobalFraction_;

            //- Worker dynamic list of donor/acceptor pairs for efficiency
            mutable donorAcceptorDynamicList* cumulativeDonorAcceptorsPtr_;

            //- Switch to control cacheing
            Switch cacheFringe_;

            //- Iteration counter
            mutable label fringeIter_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            overlapFringe(const overlapFringe&);

            //- Disallow default bitwise assignment
            void operator=(const overlapFringe&);


        // Addressing calculation functions

            //- Calculate hole-acceptor addressing
            void calcAddressing() const;

            //- Clear hole-acceptor addressing
            void clearAddressing() const;


public:

    //- Runtime type information
    TypeName("overlap");


    // Static Member Functions

        //- Helper function used to evaluate coupled boundaries, excluding
        //  overset patch
        static void evaluateNonOversetBoundaries
        (
            volScalarField::GeometricBoundaryField& psib
        );


    // Constructors

        //- Construct from dictionary
        overlapFringe
        (
            const fvMesh& mesh,
            const oversetRegion& region,
            const dictionary& dict
        );


    //- Destructor
    virtual ~overlapFringe();


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
