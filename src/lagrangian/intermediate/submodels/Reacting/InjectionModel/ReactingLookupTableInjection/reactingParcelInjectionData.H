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
    Foam::reactingParcelInjectionData

Description
    Container class to provide injection data for reacting parcels

SourceFiles
    reactingParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef reactingParcelInjectionData_H
#define reactingParcelInjectionData_H

#include "thermoParcelInjectionData.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class reactingParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const reactingParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    reactingParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
               Class reactingParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class reactingParcelInjectionData
:
    public thermoParcelInjectionData
{
protected:

    // Parcel properties

        //- List of mass fractions
        scalarList Y_;


public:

    //- Runtime type information
    TypeName("reactingParcelInjectionData");

    // Constructors

        //- Null constructor
        reactingParcelInjectionData();

        //- Construct from dictionary
        reactingParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        reactingParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~reactingParcelInjectionData();


    // Access

        //- Return const access to the list of mass fractions
        inline const scalarList& Y() const;


    // Edit

        //- Return access to the mass fractions
        inline scalarList& Y();


    // I-O

        // Friend Operators

            friend bool operator==
            (
                const reactingParcelInjectionData& a,
                const reactingParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const reactingParcelInjectionData&, "
                        "const reactingParcelInjectionData&"
                    ")"
                );

                return false;
            }

            friend bool operator!=
            (
                const reactingParcelInjectionData& a,
                const reactingParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const reactingParcelInjectionData&, "
                        "const reactingParcelInjectionData&"
                    ")"
                );

                return false;
            }

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const reactingParcelInjectionData& data
            );

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                reactingParcelInjectionData& data
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactingParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
