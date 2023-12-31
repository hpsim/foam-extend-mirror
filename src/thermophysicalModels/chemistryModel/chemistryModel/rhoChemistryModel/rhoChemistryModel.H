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
    Foam::rhoChemistryModel

Description
    Chemistry model for density-based thermodynamics

SourceFiles
    rhoChemistryModelI.H
    rhoChemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef rhoChemistryModel_H
#define rhoChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hsReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class rhoChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class rhoChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoChemistryModel(const rhoChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoChemistryModel&);


protected:

    // Protected data

        //- Thermo package
        autoPtr<hsReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("rhoChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh,
            const objectRegistry& obj,
            const word& compTypeName,
            const word& thermoTypeName
        ),
        (mesh, obj, compTypeName, thermoTypeName)
    );


    // Constructors

        //- Construct from mesh
        rhoChemistryModel
        (
            const fvMesh& mesh,
            const objectRegistry& obj,
            const word& thermoTypeName
        );


    //- Selector
    static autoPtr<rhoChemistryModel> New(const fvMesh&, const objectRegistry&);

    static autoPtr<rhoChemistryModel> New(const fvMesh& mesh)
    {
        return New(mesh, mesh);
    }


    //- Destructor
    virtual ~rhoChemistryModel();


    // Member Functions

        //- Return access to the thermo package
        inline hsReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const hsReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
