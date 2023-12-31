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
    Foam::LESmodels::SpalartAllmarasIDDES

Description
    SpalartAllmarasIDDES LES turbulence model for incompressible flows

SourceFiles
    SpalartAllmarasIDDES.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasIDDES_H
#define SpalartAllmarasIDDES_H

#include "SpalartAllmaras.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                  Class SpalartAllmarasIDDES Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmarasIDDES
:
    public SpalartAllmaras
{
    // Private data

        // Model constants

            dimensionedScalar fwStar_;
            dimensionedScalar cl_;
            dimensionedScalar ct_;


    // Private member functions

        tmp<volScalarField> alpha() const;
        tmp<volScalarField> ft(const volScalarField& S) const;
        tmp<volScalarField> fl(const volScalarField& S) const;

        tmp<volScalarField> rd
        (
            const volScalarField& visc,
            const volScalarField& S
        ) const;

        //- Delay function
        tmp<volScalarField> fd(const volScalarField& S) const;

        //- Disallow default bitwise copy construct
        SpalartAllmarasIDDES(const SpalartAllmarasIDDES&) = delete;

        //- Disallow default bitwise assignment
        SpalartAllmarasIDDES& operator=(const SpalartAllmarasIDDES&) = delete;


protected:

    // Protected member functions

        //- Length scale
        virtual tmp<volScalarField> dTilda(const volScalarField& S) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmarasIDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasIDDES
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~SpalartAllmarasIDDES() = default;


    // Member Functions

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
