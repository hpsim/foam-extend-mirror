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

Application
    dsmcFoam

Description
    Initialise a case for dsmcFoam by reading the initialisation dictionary
    system/dsmcInitialise

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOdictionary dsmcInitialiseDict
    (
        IOobject
        (
            "dsmcInitialiseDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Initialising dsmc for Time = " << runTime.timeName() << nl << endl;

    dsmcCloud dsmc("dsmc", mesh, dsmcInitialiseDict);

    label totalMolecules = dsmc.size();

    if (Pstream::parRun())
    {
        reduce(totalMolecules, sumOp<label>());
    }

    Info<< nl << "Total number of molecules added: " << totalMolecules
        << nl << endl;

    IOstream::defaultPrecision(15);

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing dsmcCloud."
            << nl << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
