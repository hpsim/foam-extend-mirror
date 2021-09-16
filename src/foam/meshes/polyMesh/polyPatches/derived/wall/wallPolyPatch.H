/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::wallPolyPatch

Description
    Wall patch.

    Added a switch for closed walls that undergo solid body motion only.

SourceFiles
    wallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wallPolyPatch_H
#define wallPolyPatch_H

#include "polyPatch.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class wallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class wallPolyPatch
:
    public polyPatch
{
    // Private data

        //- Adjust for solid body motion
        //  True for closed wall surfaces undergoing solid body motion only
        Switch closedSolidBodyMotion_;


public:

    //- Runtime type information
    TypeName("wall");


    // Constructors

        //- Construct from components
        wallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        wallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        wallPolyPatch
        (
            const wallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct as copy
        wallPolyPatch(const wallPolyPatch&);

        //- Construct as copy, resetting the boundary mesh
        wallPolyPatch(const wallPolyPatch&, const polyBoundaryMesh&);

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new wallPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new wallPolyPatch(*this, bm, index, newSize, newStart)
            );
        }

        // Virtual function for wall handling of derived class
        virtual bool isWall() const
        {
            return true;
        }

    //- Destructor
    virtual ~wallPolyPatch() = default;

    // Member functions

        // Access

            //- Return the solid body motion state
            //  True for closed wall surfaces undergoing solid body motion only
            bool closedSolidBodyMotion() const
            {
                return closedSolidBodyMotion_;
            }

            //- Return access to the solid body motion switch
            //  to allow adjustment
            //  True for closed wall surfaces undergoing solid body motion only
            Switch& closedSolidBodyMotion()
            {
                return closedSolidBodyMotion_;
            }

        //- Write the wallPolyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
