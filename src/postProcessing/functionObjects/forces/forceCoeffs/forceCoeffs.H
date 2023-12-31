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
    Foam::forceCoeffs

Description
    Derived from the forces function object, creates a specialisation to
    calculate lift and drag forces.

SourceFiles
    forceCoeffs.C
    IOforceCoeffs.H

\*---------------------------------------------------------------------------*/

#ifndef forceCoeffs_H
#define forceCoeffs_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class forceCoeffs Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffs
:
    public forces
{
    // Private data

        // Force coefficient geometry

            //- Lift
            vector liftDir_;

            //- Drag
            vector  dragDir_;

            //- Pitch
            vector pitchAxis_;


        // Free-stream conditions

            //- Velocity magnitude
            scalar magUInf_;


        // Reference scales

            //- Length
            scalar lRef_;

            //- Area
            scalar Aref_;


    // Private member functions

        //- Disallow default bitwise copy construct
        forceCoeffs(const forceCoeffs&);

        //- Disallow default bitwise assignment
        void operator=(const forceCoeffs&);


protected:

    //- Output file header information
    virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("forceCoeffs");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        forceCoeffs
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~forceCoeffs();


    // Member Functions

        //- Read the forces data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet()
        {}

        //- Write the forces
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
