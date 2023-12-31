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
    moveDynamicMesh

Description
    Mesh motion and topological mesh changes utility.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("checkFrequency", "int");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    // Read check frequency
    label checkFrequency = 1;
    args.optionReadIfPresent("checkFrequency", checkFrequency);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        if (isDir(runTime.path()/"VTK"))
        {
            Info << "Clear VTK directory" << endl;
            rmDir(runTime.path()/"VTK");
        }

        mesh.update();

#       include "volContinuity.H"
#       include "meshCourantNo.H"

        if (runTime.timeIndex() % checkFrequency == 0)
        {
            mesh.checkMesh(true);

            volScalarField magMeshCo
            (
                "magMeshCo",
                fvc::surfaceSum
                (
                    mag
                    (
                        mesh.phi()*
                        mesh.surfaceInterpolation::deltaCoeffs()/
                        mesh.magSf()
                    )
                )
            );

            if (runTime.outputTime())
            {
                Info<< "Writing mesh motion Co number" << endl;
                magMeshCo.write();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
