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
    Foam::phaseProperties

Description
    Helper class to manage multi-component phase properties

SourceFiles
    phaseProperties.C
    phasePropertiesIO.C

\*---------------------------------------------------------------------------*/

#ifndef phaseProperties_H
#define phaseProperties_H

#include "NamedEnum.H"
#include "Tuple2.H"
#include "PtrList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class phaseProperties;

Istream& operator>>(Istream&, phaseProperties&);
Ostream& operator<<(Ostream&, const phaseProperties&);


/*---------------------------------------------------------------------------*\
                      Class phaseProperties Declaration
\*---------------------------------------------------------------------------*/

class phaseProperties
{
public:

    // Public data

        //- Phase type enumeration
        enum phaseType
        {
            GAS,
            LIQUID,
            SOLID,
            UNKNOWN
        };

        //- Corresponding word representations for phase type enumerations
        static const NamedEnum<phaseType, 4> phaseTypeNames_;


private:

   // Private data

        //- Phase type
        phaseType phase_;

        //- State label (s), (l), (g) etc.
        word stateLabel_;

        //- List of component names
        wordList names_;

        //- List of component mass fractions
        scalarField Y_;

        //- Global ids
        labelList globalIds_;

        //- Map to carrier global id
        labelList globalCarrierIds_;


    // Private member functions

        //- Set global ids
        void setGlobalIds(const wordList& globalNames);

        //- Set global carrier ids - attempts to map component names to global
        //  carrier species
        void setGlobalCarrierIds(const wordList& carrierNames);

        //- Check the total mass fraction
        void checkTotalMassFraction() const;

        //- Set the state label
        word phaseToStateLabel(const phaseType pt) const;


public:

    // Constructors

        //- Null constructor
        phaseProperties();

        //- Construct from Istream
        phaseProperties(Istream&);


    //- Destructor
    ~phaseProperties() = default;


    // Public member functions

        //- Initialise the global ids
        void initialiseGlobalIds
        (
            const wordList& gasNames,
            const wordList& liquidNames,
            const wordList& solidNames
        );


        // Access

            //- Return const access to the phase type
            phaseType phase() const;

            //- Return const access to the phase state label
            const word& stateLabel() const;

            //- Return word representation of the phase type
            word phaseTypeName() const;

            //- Return the list of component names
            const wordList& names() const;

            //- Return const access to a component name
            const word& name(const label cmptI) const;

            //- Return const access to all component mass fractions
            const scalarField& Y() const;

            //- Return non-const access to a component mass fraction
            scalar& Y(const label cmptI);

            //- Return const access to the global ids
            const labelList& globalIds() const;

            //- Return const access to the map to the carrier global ids
            const labelList& globalCarrierIds() const;

            //- Return the global id of a component in the local list by name
            //  Returns -1 if not found
            label globalId(const word& cmptName) const;

            //- Return the id of a component in the local list by name
            //  Returns -1 if not found
            label id(const word& cmptName) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, phaseProperties&);
        friend Ostream& operator<<(Ostream&, const phaseProperties&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
