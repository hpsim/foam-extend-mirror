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
    Foam::XiModels::transport

Description
    Simple transport model for Xi based on Gulders correlation
    with a linear correction function to give a plausible profile for Xi.
    See report TR/HGW/10 for details on the Weller two equations model.
    See \link XiModel.H \endlink for more details on flame wrinkling modelling.

SourceFiles
    transport.C

\*---------------------------------------------------------------------------*/

#ifndef transport_H
#define transport_H

#include "XiModel.H"
#include "XiEqModel.H"
#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                          Class transport Declaration
\*---------------------------------------------------------------------------*/

class transport
:
    public XiModel
{
    // Private data

        scalar XiShapeCoef;

        autoPtr<XiEqModel> XiEqModel_;
        autoPtr<XiGModel> XiGModel_;


    // Private Member Functions

        //- Disallow copy construct
        transport(const transport&);

        //- Disallow default bitwise assignment
        void operator=(const transport&);


public:

    //- Runtime type information
    TypeName("transport");


    // Constructors

        //- Construct from components
        transport
        (
            const dictionary& XiProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~transport();


    // Member Functions

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const;

        //- Add Xi to the multivariateSurfaceInterpolationScheme table
        virtual void addXi
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
        )
        {
            fields.add(Xi_);
        }

        //- Correct the flame-wrinking Xi
        virtual void correct()
        {
            notImplemented("transport::correct()");
        }

        //- Correct the flame-wrinking Xi using the given convection scheme
        virtual void correct(const fv::convectionScheme<scalar>& mvConvection);

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
