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
    reconstructPar

Description
    Reconstructs a mesh and fields of a case that is decomposed for parallel
    execution of FOAM.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorMeshes.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "tetPointFieldReconstructor.H"
#include "reconstructLagrangian.H"

#include "faCFD.H"
#include "faMesh.H"
#include "processorFaMeshes.H"
#include "faFieldReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // enable -constant ... if someone really wants it
    // enable -noZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);
    argList::noParallel();
#   include "addRegionOption.H"
    argList::validOptions.insert("fields", "\"(list of fields)\"");
    argList::validOptions.insert("noLagrangian", "");

#   include "setRootCase.H"
#   include "createTime.H"

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    bool noLagrangian = args.optionFound("noLagrangian");

    // Determine the processor count directly
    label nProcs = 0;
    while (isDir(args.path()/(word("processor") + name(nProcs))))
    {
        ++nProcs;
    }

    if (!nProcs)
    {
        FatalErrorIn(args.executable())
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );
    }

    // use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    if (timeDirs.empty())
    {
        FatalErrorIn(args.executable())
            << "No times selected"
            << exit(FatalError);
    }

#   include "createNamedMesh.H"
    fileName regionPrefix = "";
    if (regionName != fvMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    // Set all times on processor meshes equal to reconstructed mesh
    forAll (databases, procI)
    {
        databases[procI].setTime(runTime, runTime.timeIndex());
    }

    // Read all meshes and addressing to reconstructed mesh
    processorMeshes procMeshes(databases, regionName);

    // Loop over all times
    forAll (timeDirs, timeI)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timeI], timeI);

        Info << "Time = " << runTime.timeName() << endl << endl;

        // Set time for all databases
        forAll (databases, procI)
        {
            databases[procI].setTime(timeDirs[timeI], timeI);
        }

        // Check if any new meshes need to be read.
        fvMesh::readUpdateState meshStat = mesh.readUpdate();

        fvMesh::readUpdateState procStat = procMeshes.readUpdate();

        if (procStat == fvMesh::POINTS_MOVED)
        {
            // Reconstruct the points for moving mesh cases and write them out
            procMeshes.reconstructPoints(mesh);
        }
        else if (meshStat != procStat)
        {
            WarningIn(args.executable())
                << "readUpdate for the reconstructed mesh:" << meshStat << nl
                << "readUpdate for the processor meshes  :" << procStat << nl
                << "These should be equal or your addressing"
                << " might be incorrect."
                << " Please check your time directories for any "
                << "mesh directories." << endl;
        }


        // Get list of objects from processor0 database
        IOobjectList objects(procMeshes.meshes()[0], databases[0].timeName());

        // If there are any FV fields, reconstruct them

        if
        (
            objects.lookupClass(volScalarField::typeName).size()
         || objects.lookupClass(volVectorField::typeName).size()
         || objects.lookupClass(volSphericalTensorField::typeName).size()
         || objects.lookupClass(volSymmTensorField::typeName).size()
         || objects.lookupClass(volTensorField::typeName).size()
         || objects.lookupClass(surfaceScalarField::typeName).size()
         || objects.lookupClass(surfaceVectorField::typeName).size()
         || objects.lookupClass(surfaceSphericalTensorField::typeName).size()
         || objects.lookupClass(surfaceSymmTensorField::typeName).size()
         || objects.lookupClass(surfaceTensorField::typeName).size()
        )
        {
            Info << "Reconstructing FV fields" << nl << endl;

            fvFieldReconstructor fvReconstructor
            (
                mesh,
                procMeshes.meshes(),
                procMeshes.faceProcAddressing(),
                procMeshes.cellProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            fvReconstructor.reconstructFvVolumeFields<scalar>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<vector>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<symmTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<tensor>
            (
                objects,
                selectedFields
            );

            fvReconstructor.reconstructFvSurfaceFields<scalar>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<vector>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<symmTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<tensor>
            (
                objects,
                selectedFields
            );
        }
        else
        {
            Info << "No FV fields" << nl << endl;
        }


        // If there are any point fields, reconstruct them
        if
        (
            objects.lookupClass(pointScalarField::typeName).size()
         || objects.lookupClass(pointVectorField::typeName).size()
         || objects.lookupClass(pointSphericalTensorField::typeName).size()
         || objects.lookupClass(pointSymmTensorField::typeName).size()
         || objects.lookupClass(pointTensorField::typeName).size()
        )
        {
            Info << "Reconstructing point fields" << nl << endl;

            pointMesh pMesh(mesh);
            PtrList<pointMesh> pMeshes(procMeshes.meshes().size());

            forAll (pMeshes, procI)
            {
                pMeshes.set(procI, new pointMesh(procMeshes.meshes()[procI]));
            }

            pointFieldReconstructor pointReconstructor
            (
                pMesh,
                pMeshes,
                procMeshes.pointProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            pointReconstructor.reconstructFields<scalar>(objects);
            pointReconstructor.reconstructFields<vector>(objects);
            pointReconstructor.reconstructFields<sphericalTensor>(objects);
            pointReconstructor.reconstructFields<symmTensor>(objects);
            pointReconstructor.reconstructFields<tensor>(objects);
        }
        else
        {
            Info << "No point fields" << nl << endl;
        }

        // If there are any tetFem fields, reconstruct them
        if
        (
            objects.lookupClass(tetPointScalarField::typeName).size()
         || objects.lookupClass(tetPointVectorField::typeName).size()
         || objects.lookupClass(tetPointSphericalTensorField::typeName).size()
         || objects.lookupClass(tetPointSymmTensorField::typeName).size()
         || objects.lookupClass(tetPointTensorField::typeName).size()

         || objects.lookupClass(elementScalarField::typeName).size()
         || objects.lookupClass(elementVectorField::typeName).size()
        )
        {
            Info << "Reconstructing tet point fields" << nl << endl;

            tetPolyMesh tetMesh(mesh);
            PtrList<tetPolyMesh> tetMeshes(procMeshes.meshes().size());

            forAll (tetMeshes, procI)
            {
                tetMeshes.set
                (
                    procI,
                    new tetPolyMesh(procMeshes.meshes()[procI])
                );
            }

            tetPointFieldReconstructor tetPointReconstructor
            (
                tetMesh,
                tetMeshes,
                procMeshes.pointProcAddressing(),
                procMeshes.faceProcAddressing(),
                procMeshes.cellProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            // Reconstruct tet point fields
            tetPointReconstructor.reconstructTetPointFields<scalar>(objects);
            tetPointReconstructor.reconstructTetPointFields<vector>(objects);
            tetPointReconstructor.
                reconstructTetPointFields<sphericalTensor>(objects);
            tetPointReconstructor.
                reconstructTetPointFields<symmTensor>(objects);
            tetPointReconstructor.reconstructTetPointFields<tensor>(objects);

            tetPointReconstructor.reconstructElementFields<scalar>(objects);
            tetPointReconstructor.reconstructElementFields<vector>(objects);
        }
        else
        {
            Info << "No tetFem fields" << nl << endl;
        }


        // If there are any clouds, reconstruct them.
        // The problem is that a cloud of size zero will not get written so
        // in pass 1 we determine the cloud names and per cloud name the
        // fields. Note that the fields are stored as IOobjectList from
        // the first processor that has them. They are in pass2 only used
        // for name and type (scalar, vector etc).

        if (!noLagrangian)
        {
            HashTable<IOobjectList> cloudObjects;

            forAll (databases, procI)
            {
                fileNameList cloudDirs
                (
                    readDir
                    (
                        databases[procI].timePath()/regionPrefix/cloud::prefix,
                        fileName::DIRECTORY
                    )
                );

                forAll (cloudDirs, i)
                {
                    // Check if we already have cloud objects for this cloudname
                    HashTable<IOobjectList>::const_iterator iter =
                        cloudObjects.find(cloudDirs[i]);

                    if (iter == cloudObjects.end())
                    {
                        // Do local scan for valid cloud objects
                        IOobjectList sprayObjs
                        (
                            procMeshes.meshes()[procI],
                            databases[procI].timeName(),
                            cloud::prefix/cloudDirs[i]
                        );

                        IOobject* positionsPtr = sprayObjs.lookup("positions");

                        if (positionsPtr)
                        {
                            cloudObjects.insert(cloudDirs[i], sprayObjs);
                        }
                    }
                }
            }


            if (cloudObjects.size())
            {
                // Pass2: reconstruct the cloud
                forAllConstIter(HashTable<IOobjectList>, cloudObjects, iter)
                {
                    const word cloudName = string::validate<word>(iter.key());

                    // Objects (on arbitrary processor)
                    const IOobjectList& sprayObjs = iter();

                    Info<< "Reconstructing lagrangian fields for cloud "
                        << cloudName << nl << endl;

                    reconstructLagrangianPositions
                    (
                        mesh,
                        cloudName,
                        procMeshes.meshes(),
                        procMeshes.faceProcAddressing(),
                        procMeshes.cellProcAddressing()
                    );
                    reconstructLagrangianFields<label>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<scalar>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<vector>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<sphericalTensor>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<symmTensor>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<tensor>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                }
            }
            else
            {
                Info << "No lagrangian fields" << nl << endl;
            }
        }

        // If there are any FA fields, reconstruct them

        if
        (
            objects.lookupClass(areaScalarField::typeName).size()
         || objects.lookupClass(areaVectorField::typeName).size()
         || objects.lookupClass(areaSphericalTensorField::typeName).size()
         || objects.lookupClass(areaSymmTensorField::typeName).size()
         || objects.lookupClass(areaTensorField::typeName).size()
         || objects.lookupClass(edgeScalarField::typeName).size()
        )
        {
            Info << "Reconstructing FA fields" << nl << endl;

            faMesh aMesh(mesh);

            // Create mesh addressing by reading from files
            processorFaMeshes procFaMeshes(procMeshes.meshes(), true);

            faFieldReconstructor faReconstructor
            (
                aMesh,
                procFaMeshes.meshes(),
                procFaMeshes.edgeProcAddressing(),
                procFaMeshes.faceProcAddressing(),
                procFaMeshes.boundaryProcAddressing()
            );

            faReconstructor.reconstructFaAreaFields<scalar>(objects);
            faReconstructor.reconstructFaAreaFields<vector>(objects);
            faReconstructor
               .reconstructFaAreaFields<sphericalTensor>(objects);
            faReconstructor.reconstructFaAreaFields<symmTensor>(objects);
            faReconstructor.reconstructFaAreaFields<tensor>(objects);

            faReconstructor.reconstructFaEdgeFields<scalar>(objects);
        }
        else
        {
            Info << "No FA fields" << nl << endl;
        }

        // If there are any "uniform" directories copy them from
        // the master processor

        fileName uniformDir0 = databases[0].timePath()/"uniform";
        if (isDir(uniformDir0))
        {
            cp(uniformDir0, runTime.timePath());
        }
    }

    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
