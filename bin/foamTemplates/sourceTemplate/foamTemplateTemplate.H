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
    Foam::ClassName

Description

Author

SourceFiles
    ClassNameI.H
    ClassName.C
    ClassNameIO.C

\*---------------------------------------------------------------------------*/

#ifndef ClassName_H
#define ClassName_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class someClass;


// Forward declaration of friend functions and operators

template<TemplateClassArgument>
class ClassName;

template<TemplateClassArgument>
Istream& operator>>(Istream&, ClassName<TemplateArgument>&);

template<TemplateClassArgument>
Ostream& operator<<(Ostream&, const ClassName<TemplateArgument>&);


/*---------------------------------------------------------------------------*\
                         Class ClassName Declaration
\*---------------------------------------------------------------------------*/

template<TemplateClassArgument>
class ClassName
:
    public baseClassName
{
    // Private data

        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ClassName(const ClassName<TemplateArgument>&);

        //- Disallow default bitwise assignment
        void operator=(const ClassName<TemplateArgument>&);


public:

    // Static data members

        //- Static data someStaticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        ClassName();

        //- Construct from components
        ClassName(const dataType& data);

        //- Construct from Istream
        ClassName(Istream&);

        //- Construct as copy
        ClassName(const ClassName<TemplateArgument>&);


    // Selectors

        //- Select null constructed
        static autoPtr<ClassName<TemplateArgument> > New();


    //- Destructor
        ~ClassName();


    // Member Functions

    // Member Operators

        void operator=(const ClassName<TemplateArgument>&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>> <TemplateArgument>
        (Istream&, ClassName<TemplateArgument>&);

        friend Ostream& operator<< <TemplateArgument>
        (Ostream&, const ClassName<TemplateArgument>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ClassNameI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ClassName.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
