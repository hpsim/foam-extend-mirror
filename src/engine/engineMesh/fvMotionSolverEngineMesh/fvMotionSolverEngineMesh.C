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

\*---------------------------------------------------------------------------*/

#include "fvMotionSolverEngineMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fvMotionSolverEngineMesh, 0);

addToRunTimeSelectionTable(engineMesh, fvMotionSolverEngineMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvMotionSolverEngineMesh::fvMotionSolverEngineMesh(const IOobject& io)
:
    engineMesh(io),
    pistonLayers_("pistonLayers", dimLength, 0.0),
    motionSolver_(*this, engineDB_.engineDict().lookup("motionSolver"))
{
    if (engineDB_.engineDict().found("pistonLayers"))
    {
        engineDB_.engineDict().lookup("pistonLayers") >> pistonLayers_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fvMotionSolverEngineMesh::~fvMotionSolverEngineMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fvMotionSolverEngineMesh::move()
{
    scalar deltaZ = engineDB_.pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << endl;

    // Position of the top of the static mesh layers above the piston
    scalar pistonPlusLayers = pistonPosition_.value() + pistonLayers_.value();

    motionSolver_.cellMotionU().boundaryField()[pistonIndex_] == deltaZ;

    {
        scalarField linerPoints =
            motionSolver_.cellMotionU()
           .boundaryField()[linerIndex_].patch().Cf().component(vector::Z);

        motionSolver_.cellMotionU().boundaryField()[linerIndex_] ==
            deltaZ*pos(deckHeight_.value() - linerPoints)
            *(deckHeight_.value() - linerPoints)
            /(deckHeight_.value() - pistonPlusLayers);
    }

    motionSolver_.solve();

    if (engineDB_.foundObject<surfaceScalarField>("phi"))
    {
        surfaceScalarField& phi =
            const_cast<surfaceScalarField&>
            (engineDB_.lookupObject<surfaceScalarField>("phi"));

        const volScalarField& rho =
            engineDB_.lookupObject<volScalarField>("rho");

        const volVectorField& U =
            engineDB_.lookupObject<volVectorField>("U");

        bool absolutePhi = false;
        if (moving())
        {
            phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);
            absolutePhi = true;
        }

        movePoints(motionSolver_.curPoints());

        if (absolutePhi)
        {
            phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);
        }
    }
    else
    {
        movePoints(motionSolver_.curPoints());
    }


    pistonPosition_.value() += deltaZ;
    scalar pistonSpeed = deltaZ/engineDB_.deltaT().value();

    Info<< "clearance: " << deckHeight_.value() - pistonPosition_.value() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
