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
    Foam::sampledIsoSurfaceCell

Description
    A sampledSurface defined by a surface of iso value. Always triangulated.
    To be used in sampleSurfaces / functionObjects. Recalculates iso surface
    only if time changes.

SourceFiles
    sampledIsoSurfaceCell.C

\*---------------------------------------------------------------------------*/

#ifndef sampledIsoSurfaceCell_H
#define sampledIsoSurfaceCell_H

#include "sampledSurface.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class sampledIsoSurfaceCell Declaration
\*---------------------------------------------------------------------------*/

class sampledIsoSurfaceCell
:
    public sampledSurface,
    public triSurface
{
    // Private data

        //- Field to get isoSurface of
        const word isoField_;

        //- iso value
        const scalar isoVal_;

        //- Whether to coarse
        const Switch regularise_;

        //- Whether to recalculate cell values as average of point values
        const Switch average_;

        //- zone name (if restricted to zones)
        word zoneName_;

        //- triangles converted to faceList
        mutable autoPtr<faceList> facesPtr_;


        // Recreated for every isoSurface

            //- Time at last call, also track it surface needs an update
            mutable label prevTimeIndex_;

            //- For every triangle the original cell in mesh
            mutable labelList meshCells_;


    // Private Member Functions

        //- Create iso surface (if time has changed)
        //  Do nothing (and return false) if no update was needed
        bool updateGeometry() const;

        //- sample field on faces
        template <class Type>
        tmp<Field<Type> > sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        template <class Type>
        tmp<Field<Type> >
        interpolateField(const interpolation<Type>&) const;


public:

    //- Runtime type information
    TypeName("sampledIsoSurfaceCell");


    // Constructors

        //- Construct from dictionary
        sampledIsoSurfaceCell
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~sampledIsoSurfaceCell();


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();


        //- Points of surface
        virtual const pointField& points() const
        {
            return triSurface::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            if (facesPtr_.empty())
            {
                const triSurface& s = *this;

                facesPtr_.reset(new faceList(s.size()));

                forAll(s, i)
                {
                    facesPtr_()[i] = s[i].triFaceFace();
                }
            }
            return facesPtr_;
        }


        //- sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const;

        //- sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const;

        //- sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const;

        //- sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const;

        //- sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const;


        //- interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const;

        //- interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const;

        //- interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const;

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sampledIsoSurfaceCellTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
