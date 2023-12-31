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
    Foam::heatFlux

Description

SourceFiles
    heatFlux.C
    heatFluxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef heatFlux_H
#define heatFlux_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"
#include "Switch.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class heatFluxItem;
class OFstream;
template<class Type>
class List;

/*---------------------------------------------------------------------------*\
                         Class heatFlux Declaration
\*---------------------------------------------------------------------------*/

class heatFlux
{
    // Private data

        //- Name of this set of heat flux check
        word name_;

        //- name of fluid thermal diffusivity
        word Kfluid_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;


    // Private Member Functions

        // Calculation functions

            //- Write
            virtual void calcAndPrint();


        //- Disallow default bitwise copy construct
        heatFlux(const heatFlux&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const heatFlux&) = delete;


public:

    //- Runtime type information
    TypeName("heatFlux");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        heatFlux
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~heatFlux() = default;


    // Member Functions

        //- Return name of the set of field averages
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field average data
        virtual void read(const dictionary&);

        //- Execute the averaging
        virtual void execute();

        //- Execute the averaging at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the field average data and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const pointField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
