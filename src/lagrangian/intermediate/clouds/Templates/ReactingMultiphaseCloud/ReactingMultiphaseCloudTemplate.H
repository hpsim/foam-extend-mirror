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
    Foam::ReactingMultiphaseCloud

Description
    Templated base class for multiphase reacting cloud

    - Adds to reacting cloud
      - multiphase composition
      - devolatilisatsion
      - surface reactions

SourceFiles
    ReactingMultiphaseCloudTemplateI.H
    ReactingMultiphaseCloudTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseCloudTemplate_H
#define ReactingMultiphaseCloudTemplate_H

#include "ReactingCloudTemplate.H"
#include "reactingMultiphaseCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DevolatilisationModel;

template<class CloudType>
class SurfaceReactionModel;

/*---------------------------------------------------------------------------*\
                   Class ReactingMultiphaseCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingMultiphaseCloud
:
    public ReactingCloud<ParcelType>,
    public reactingMultiphaseCloud
{
    // Private member functions

        //- Disallow default bitwise copy construct
        ReactingMultiphaseCloud(const ReactingMultiphaseCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ReactingMultiphaseCloud&);


protected:

    // Protected data

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Devolatilisation model
            autoPtr
            <
                DevolatilisationModel<ReactingMultiphaseCloud<ParcelType> >
            >
            devolatilisationModel_;

            //- Surface reaction model
            autoPtr
            <
                SurfaceReactionModel<ReactingMultiphaseCloud<ParcelType> >
            >
            surfaceReactionModel_;


        // Check

            //- Total mass transferred to continuous phase via devolatilisation
            scalar dMassDevolatilisation_;

            //- Total mass transferred to continuous phase via surface
            //  reactions
            scalar dMassSurfaceReaction_;


    // Protected member functions

        // Cloud evolution functions

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolveCloud();

            //- Post-evolve
            void postEvolve();


public:

    // Constructors

        //- Construct given carrier gas fields
        ReactingMultiphaseCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            basicThermo& thermo,
            bool readFields = true
        );


    //- Destructor
    virtual ~ReactingMultiphaseCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            //- Return the constant properties
            inline const typename ParcelType::constantProperties&
                constProps() const;


            // Sub-models

                //- Return reference to devolatilisation model
                inline const DevolatilisationModel
                <
                    ReactingMultiphaseCloud<ParcelType>
                >&
                devolatilisation() const;

                //- Return reference to reacting surface reaction model
                inline const SurfaceReactionModel
                <
                    ReactingMultiphaseCloud<ParcelType>
                >&
                surfaceReaction() const;


        // Check

            //- Print cloud information
            void info() const;

            //- Add to cumulative volatilisation mass transfer
            void addToMassDevolatilisation(const scalar dMass);

            //- Add to cumulative surface reaction transfer
            void addToMassSurfaceReaction(const scalar dMass);


        // Cloud evolution functions

            //- Check parcel properties
            void checkParcelProperties
            (
                ParcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (inject, move)
            void evolve();

        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseCloudTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingMultiphaseCloudTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
