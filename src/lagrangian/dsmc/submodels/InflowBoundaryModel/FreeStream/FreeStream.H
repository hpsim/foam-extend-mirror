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
    Foam::FreeStream

Description
    Inserting new particles across the faces of a all patched of type
    "patch" for a free stream.  Uniform values number density, temperature
    and velocity sourced face-by-face from the boundaryT and boundaryU fields
    of the cloud.

\*---------------------------------------------------------------------------*/

#ifndef FreeStream_H
#define FreeStream_H

#include "InflowBoundaryModel.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class FreeStream Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FreeStream
:
    public InflowBoundaryModel<CloudType>
{
    // Private data

        //- The indices of patches to introduce molecules across
        labelList patches_;

        //- The molecule types to be introduced
        labelList moleculeTypeIds_;

        //- The number density of the species in the inflow
        Field<scalar> numberDensities_;

        //- A List of Lists of Fields specifying carry-over of mass flux from
        // one timestep to the next
        // + Outer List - one inner List for each patch
        // + Inner List - one Field for every species to be introduced
        // + Each field entry corresponding to a face to be injected across
        //   with a particular species
        List<List<Field<scalar> > > particleFluxAccumulators_;


public:

    //- Runtime type information
    TypeName("FreeStream");


    // Constructors

        //- Construct from dictionary
        FreeStream
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Destructor
    virtual ~FreeStream();


    // Member Functions

        //- Introduce particles
        virtual void inflow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FreeStream.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
