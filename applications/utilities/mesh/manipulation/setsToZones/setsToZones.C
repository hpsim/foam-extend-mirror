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

Description
    Add pointZones/faceZones/cellZones to the mesh from similar named
    pointSets/faceSets/cellSets.

    For faceZones the user needs to specify a flip
    condition which denotes the side of the face.  This application
    reads a cellSet (xxxMasterCells if 'xxx' is the name of the faceSet) which
    is the masterCells of the zone.  Master cell is the one IN FRONT of the
    face, ie. the one into which the face normal points.  If master cells are
    not found, take faces without a flip

    If one is not interested in sidedness specify the -noFlipMap
    command line option.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "IStringStream.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IOobjectList.H"
#include "SortableList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("noFlipMap", "");

#   include "addRegionOption.H"
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    bool noFlipMap = args.optionFound("noFlipMap");

    // Get times list
    instantList Times = runTime.times();

    label startTime = Times.size()-1;
    label endTime = Times.size();

    // check -time and -latestTime options
#   include "checkTimeOption.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createNamedPolyMesh.H"

    // Search for list of objects for the time of the mesh
    IOobjectList objects
    (
        mesh,
        mesh.facesInstance(),
        polyMesh::meshSubDir/"sets"
    );

    Info<< "Searched : " << mesh.facesInstance()/polyMesh::meshSubDir/"sets"
        << nl
        << "Found    : " << objects.names() << nl
        << endl;


    IOobjectList pointObjects(objects.lookupClass(pointSet::typeName));

    for
    (
        IOobjectList::const_iterator iter = pointObjects.begin();
        iter != pointObjects.end();
        ++iter
    )
    {
        // Set not in memory. Load it.
        pointSet set(*iter());
        SortableList<label> pointLabels(set.toc());

        label zoneID = mesh.pointZones().findZoneID(set.name());
        if (zoneID == -1)
        {
            Info<< "Adding set " << set.name() << " as a pointZone." << endl;

            label sz = mesh.pointZones().size();
            mesh.pointZones().setSize(sz+1);
            mesh.pointZones().set
            (
                sz,
                new pointZone
                (
                    set.name(),             //name
                    pointLabels,            //addressing
                    sz,                     //index
                    mesh.pointZones()       //pointZoneMesh
                )
            );
            mesh.pointZones().writeOpt() = IOobject::AUTO_WRITE;
            mesh.pointZones().instance() = mesh.facesInstance();
        }
        else
        {
            Info<< "Overwriting contents of existing pointZone " << zoneID
                << " with that of set " << set.name() << "." << endl;
            mesh.pointZones()[zoneID] = pointLabels;
            mesh.pointZones().writeOpt() = IOobject::AUTO_WRITE;
            mesh.pointZones().instance() = mesh.facesInstance();
        }
    }



    IOobjectList faceObjects(objects.lookupClass(faceSet::typeName));

    HashSet<word> masterCellSets;

    for
    (
        IOobjectList::const_iterator iter = faceObjects.begin();
        iter != faceObjects.end();
        ++iter
    )
    {
        // Set not in memory. Load it.
        faceSet set(*iter());
        SortableList<label> faceLabels(set.toc());

        dynamicLabelList addressing(set.size());
        DynamicList<bool> flipMap(set.size());

        if (!noFlipMap)
        {
            word setName(set.name() + "MasterCells");

            Info<< "Using cellSet " << setName
                << " to determine the master side of the face zone "
                << set.name() << nl
                << endl;

            // Load corresponding cells
            cellSet cells
            (
                mesh,
                setName,
                IOobject::READ_IF_PRESENT
            );

            // Store setName to exclude from cellZones further on
            masterCellSets.insert(setName);

            if (!cells.empty())
            {
                forAll(faceLabels, i)
                {
                    label faceI = faceLabels[i];

                    bool flip = false;

                    if (mesh.isInternalFace(faceI))
                    {
                        if
                        (
                            cells.found(mesh.faceOwner()[faceI])
                        && !cells.found(mesh.faceNeighbour()[faceI])
                        )
                        {
                            // Fixed, using master zone.  HJ, 17/Feb/2011
                            flip = true;
                        }
                        else if
                        (
                           !cells.found(mesh.faceOwner()[faceI])
                         && cells.found(mesh.faceNeighbour()[faceI])
                        )
                        {
                            // Fixed, using master zone.  HJ, 17/Feb/2011
                            flip = false;
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "Owner or neighbour of internal face "
                                << faceI << " should be in cellSet "
                                << cells.name()
                                << " to be able to determine orientation."
                                << nl << "Face:" << faceI
                                << " own:" << mesh.faceOwner()[faceI]
                                << " OwnInCellSet:"
                                << cells.found(mesh.faceOwner()[faceI])
                                << " nei:" << mesh.faceNeighbour()[faceI]
                                << " NeiInCellSet:"
                                << cells.found(mesh.faceNeighbour()[faceI])
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        if (cells.found(mesh.faceOwner()[faceI]))
                        {
                            // Fixed, using master zone.  HJ, 17/Feb/2011
                            flip = true;
                        }
                        else
                        {
                            // Fixed, using master zone.  HJ, 17/Feb/2011
                            flip = false;
                        }
                    }

                    addressing.append(faceI);
                    flipMap.append(flip);
                }
            }
            else
            {
                // Cell set not found or empty.  Using faces without flip
                Info<< "cellSet " << setName
                    << " not found or empty.  Setting flipMap to false" << nl
                    << endl;

                forAll(faceLabels, i)
                {
                    label faceI = faceLabels[i];
                    addressing.append(faceI);
                    flipMap.append(false);
                }
            }
        }
        else
        {
            // No flip map.
            forAll(faceLabels, i)
            {
                label faceI = faceLabels[i];
                addressing.append(faceI);
                flipMap.append(false);
            }
        }

        label zoneID = mesh.faceZones().findZoneID(set.name());
        if (zoneID == -1)
        {
            Info<< "Adding set " << set.name() << " as a faceZone." << endl;
            label sz = mesh.faceZones().size();
            mesh.faceZones().setSize(sz+1);
            mesh.faceZones().set
            (
                sz,
                new faceZone
                (
                    set.name(),             // name
                    addressing.shrink(),    // addressing
                    flipMap.shrink(),       // flipmap
                    sz,                     // index
                    mesh.faceZones()        // faceZoneMesh
                )
            );
            mesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
            mesh.faceZones().instance() = mesh.facesInstance();
        }
        else
        {
            Info<< "Overwriting contents of existing faceZone " << zoneID
                << " with that of set " << set.name() << "." << endl;
            mesh.faceZones()[zoneID].resetAddressing
            (
                addressing.shrink(),
                flipMap.shrink()
            );
            mesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
            mesh.faceZones().instance() = mesh.facesInstance();
        }
    }



    IOobjectList cellObjects(objects.lookupClass(cellSet::typeName));

    for
    (
        IOobjectList::const_iterator iter = cellObjects.begin();
        iter != cellObjects.end();
        ++iter
    )
    {
        if (!masterCellSets.found(iter.key()))
        {
            // Set not in memory. Load it.
            cellSet set(*iter());
            SortableList<label> cellLabels(set.toc());

            label zoneID = mesh.cellZones().findZoneID(set.name());
            if (zoneID == -1)
            {
                Info<< "Adding set " << set.name() << " as a cellZone."
                    << endl;

                label sz = mesh.cellZones().size();
                mesh.cellZones().setSize(sz+1);
                mesh.cellZones().set
                (
                    sz,
                    new cellZone
                    (
                        set.name(),             // name
                        cellLabels,             // addressing
                        sz,                     // index
                        mesh.cellZones()        // pointZoneMesh
                    )
                );
                mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
                mesh.cellZones().instance() = mesh.facesInstance();
            }
            else
            {
                Info<< "Overwriting contents of existing cellZone " << zoneID
                    << " with that of set " << set.name() << "." << endl;
                mesh.cellZones()[zoneID] = cellLabels;
                mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
                mesh.cellZones().instance() = mesh.facesInstance();
            }
        }
    }


    Info<< "Writing mesh." << endl;

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
