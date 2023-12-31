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
    Foam::radiation::fvDOM

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media, not including scatter.

    Available absorption models:
        greyMeanAbsoprtionEmission
        wideBandAbsorptionEmission

    i.e. dictionary
    fvDOMCoeffs
    {
        nPhi    1;          // azimuthal angles in PI/2 on X-Y.(from Y to X)
        nTheta  2;          // polar angles in PI (from Z to X-Y plane)
        convergence 1e-4;   // convergence criteria for radiation iteration
    }

    solverFreq   1; // Number of flow iterations per radiation iteration

    The total number of solid angles is  4*nPhi*nTheta.

    In 1D the direction of the rays is X (nPhi and nTheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)

SourceFiles
    fvDOM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelfvDOM_H
#define radiationModelfvDOM_H

#include "radiativeIntensityRay.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class fvDOM Declaration
\*---------------------------------------------------------------------------*/

class fvDOM
:
    public radiationModel
{
    // Private data

        //- Incident radiation  [W/m2]
        volScalarField G_;

        //- Total radiative heat flux [W/m2]
        volScalarField Qr_;

        //- Total absorption coefficient [1/m]
        volScalarField a_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nLambda_;

        //- Wavelength total absorption coefficient [1/m]
        PtrList<volScalarField> aLambda_;

        //- Black body
        blackBodyEmission blackBody_;

        //- List of pointers to radiative intensity rays
        PtrList<radiativeIntensityRay> IRay_;

        //- Emmited radiative heat flux [W/m2], per band
        PtrList<volScalarField::GeometricBoundaryField> Qem_;

        //- Incidet radiative heat flux [W/m2], per band
        PtrList<volScalarField::GeometricBoundaryField> Qin_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- List of cached fvMatrices for rays
        List<PtrList<fvScalarMatrix> >fvRayDiv_;

        //- Cache convection div matrix
        bool cacheDiv_;

        //- Maximum omega weight
        scalar omegaMax_;


    // Private member functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        fvDOM(const fvDOM&);

        //- Disallow default bitwise assignment
        void operator=(const fvDOM&);

        //- Update nlack body emission
        void updateBlackBodyEmission();


public:

    //- Runtime type information
    TypeName("fvDOM");


    // Constructors

        //- Construct from components
        fvDOM(const volScalarField& T);

        //- Construct from components
        fvDOM(const word& type, const volScalarField& T);

        //- Construct from components
        fvDOM(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~fvDOM();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            virtual bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

            //- Set the rayId and lambdaId from by decomposing an intensity
            //  field name
            void setRayIdLambdaId
            (
                const word& name,
                label& rayId,
                label& lambdaId
            ) const;

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;


        // Access

            //- Ray intensity for rayI
            inline const radiativeIntensityRay& IRay(const label rayI) const;

            //- Ray intensity for rayI and lambda bandwidth
            inline const volScalarField& IRayLambda
            (
                const label rayI,
                const label lambdaI
            ) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of wavelengths
            inline label nLambda() const;

            //- Const access to total absorption coefficient
            inline const volScalarField& a() const;

            //- Const access to wavelength total absorption coefficient
            inline const volScalarField& aLambda(const label lambdaI) const;

            //- Const access to incident radiation field
            inline const volScalarField& G() const;

            //- Const access to total radiative heat flux field
            inline const volScalarField& Qr() const;

            //- Const access to band wise emitted radiative heat flux field
            inline const volScalarField::GeometricBoundaryField& Qem
            (
                const label lambda
            ) const;

            //- Const access to band wise incident radiative heat flux field
            inline const volScalarField::GeometricBoundaryField& Qin
            (
                const label lambda
            ) const;

            //- Total incident radiative heat flux field
            tmp<volScalarField::GeometricBoundaryField> Qin() const;

            //- Total emitted radiative heat flux field
            tmp<volScalarField::GeometricBoundaryField> Qem() const;

            //- Const access to black body
            inline const blackBodyEmission& blackBody() const;

            //- Const access to cached fvMatrix
            inline const fvScalarMatrix& fvRayDiv
            (
                const label lambdaI,
                const label rayId
            ) const;

            //- Caching div(Ji, Ilamda)
            inline bool cacheDiv() const;

            //- Return omegaMax
            inline scalar omegaMax() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvDOMI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
