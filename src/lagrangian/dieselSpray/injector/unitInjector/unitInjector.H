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
    Foam::unitInjector

Description
    The unit injector

SourceFiles
    unitInjectorI.H
    unitInjector.C
    unitInjectorIO.C

\*---------------------------------------------------------------------------*/

#ifndef unitInjector_H
#define unitInjector_H

#include "injectorType.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class unitInjector Declaration
\*---------------------------------------------------------------------------*/

class unitInjector
:
    public injectorType
{

private:

    typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

    // Private data

        //- Injector properties dictionary
        dictionary propsDict_;

        //- Injector position
        vector position_;

        //- Injector direction
        vector direction_;

        //- Injector diameter
        scalar d_;

        //- Cd
        scalar Cd_;

        //- Parcel mass
        scalar mass_;

        //- Number of parcels
        label nParcels_;

        scalarField X_;
        List<pair> massFlowRateProfile_;
        List<pair> velocityProfile_;
        List<pair> injectionPressureProfile_;
        List<pair> CdProfile_;
        List<pair> TProfile_;
        scalar averageParcelMass_;

        bool pressureIndependentVelocity_;

        //- two orthogonal vectors that are also orthogonal
        //  to the injection direction
        vector tangentialInjectionVector1_, tangentialInjectionVector2_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        unitInjector(const unitInjector&);

        //- Disallow default bitwise assignment
        void operator=(const unitInjector&);

        //- Create two vectors orthonoal to each other
        //  and the injection vector
        void setTangentialVectors();

        //- Return the fraction of the total injected liquid
        scalar fractionOfInjection(const scalar time) const;


public:

    //- Runtime type information
    TypeName("unitInjector");


    // Constructors

        //- Construct from components
        unitInjector
        (
            const Time& t,
            const dictionary& dict
        );


    // Destructor

        ~unitInjector();


    // Member Functions

        //- Return number of particles to inject
        label nParcelsToInject
        (
            const scalar t0,
            const scalar t1
        ) const;

        //- Return the injection position
        const vector position(const label n) const;

        //- Return the injection position
        vector position
        (
            const label n,
            const scalar time,
            const bool twoD,
            const scalar angleOfWedge,
            const vector& axisOfSymmetry,
            const vector& axisOfWedge,
            const vector& axisOfWedgeNormal,
            Random& rndGen
        ) const;

        //- Return the number of holes
        label nHoles() const;

        //- Return the injector diameter
        scalar d() const;

        //- Return the injection direction
        const vector& direction
        (
            const label i,
            const scalar time
        ) const;

        //- Return the mass of the injected particle
        scalar mass
        (
            const scalar t0,
            const scalar t1,
            const bool twoD,
            const scalar angleOfWedge
        ) const;

        //- Return the mass injected by the injector
        scalar mass() const;

        //- Return the fuel mass fractions of the injected particle
        const scalarField& X() const;

        //- Return the temperature profile of the injected particle
        List<pair> T() const;

        //- Return the temperature of the injected particle
        scalar T(const scalar time) const;

        //- Return the start-of-injection time
        scalar tsoi() const;

        //- Return the end-of-injection time
        scalar teoi() const;

        //- Return the injected liquid mass
        scalar injectedMass(const scalar t) const;

        List<pair> massFlowRateProfile() const
        {
            return massFlowRateProfile_;
        }

        scalar massFlowRate(const scalar time) const;

        List<pair> injectionPressureProfile() const
        {
            return injectionPressureProfile_;
        }

        scalar injectionPressure(const scalar time) const;

        List<pair> velocityProfile() const
        {
            return velocityProfile_;
        }

        scalar velocity(const scalar time) const;

        List<pair> CdProfile() const;
        scalar Cd(const scalar time) const;

        vector tan1(const label n) const;
        vector tan2(const label n) const;

        void correctProfiles
        (
            const liquidMixture& fuel,
            const scalar referencePressure
        );

        bool pressureIndependentVelocity() const
        {
            return pressureIndependentVelocity_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
