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
    Tries to figure out what the refinement level is on refined cartesian
    meshes. Run BEFORE snapping.

    Writes
    - volScalarField 'refinementLevel' with current refinement level.
    - cellSet 'refCells' which are the cells that need to be refined to satisfy
      2:1 refinement.

    Works by dividing cells into volume bins.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "SortableList.H"
#include "labelIOList.H"
#include "fvMesh.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return true if any cells had to be split to keep a difference between
// neighbouring refinement levels < limitDiff. Puts cells into refCells and
// update refLevel to account for refinement.
bool limitRefinementLevel
(
    const primitiveMesh& mesh,
    labelList& refLevel,
    cellSet& refCells
)
{
    const labelListList& cellCells = mesh.cellCells();

    label oldNCells = refCells.size();

    forAll(cellCells, cellI)
    {
        const labelList& cCells = cellCells[cellI];

        forAll(cCells, i)
        {
            if (refLevel[cCells[i]] > (refLevel[cellI]+1))
            {
                // Found neighbour with >=2 difference in refLevel.
                refCells.insert(cellI);
                refLevel[cellI]++;
                break;
            }
        }
    }

    if (refCells.size() > oldNCells)
    {
        Info<< "Added an additional " << refCells.size() - oldNCells
            << " cells to satisfy 1:2 refinement level"
            << endl;

        return true;
    }
    else
    {
        return false;
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("readLevel", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    Info<< "Dividing cells into bins depending on cell volume.\nThis will"
        << " correspond to refinement levels for a mesh with only 2x2x2"
        << " refinement\n"
        << "The upper range for every bin is always 1.1 times the lower range"
        << " to allow for some truncation error."
        << nl << endl;

    bool readLevel = args.optionFound("readLevel");

    const scalarField& vols = mesh.cellVolumes();

    SortableList<scalar> sortedVols(vols);

    // All cell labels, sorted per bin.
    DynamicList<dynamicLabelList > bins;

    // Lower/upper limits
    DynamicList<scalar> lowerLimits;
    DynamicList<scalar> upperLimits;

    // Create bin0. Have upperlimit as factor times lowerlimit.
    bins.append(dynamicLabelList());
    lowerLimits.append(sortedVols[0]);
    upperLimits.append(1.1*lowerLimits[lowerLimits.size()-1]);

    forAll(sortedVols, i)
    {
        if (sortedVols[i] > upperLimits[upperLimits.size()-1])
        {
            // New value outside of current bin

            // Shrink old bin.
            dynamicLabelList& bin = bins[bins.size()-1];

            bin.shrink();

            Info<< "Collected " << bin.size() << " elements in bin "
                << lowerLimits[lowerLimits.size()-1] << " .. "
                << upperLimits[upperLimits.size()-1] << endl;

            // Create new bin.
            bins.append(dynamicLabelList());
            lowerLimits.append(sortedVols[i]);
            upperLimits.append(1.1*lowerLimits[lowerLimits.size()-1]);

            Info<< "Creating new bin " << lowerLimits[lowerLimits.size()-1]
                << " .. " << upperLimits[upperLimits.size()-1]
                << endl;
        }

        // Append to current bin.
        dynamicLabelList& bin = bins[bins.size()-1];

        bin.append(sortedVols.indices()[i]);
    }
    Info<< endl;

    bins[bins.size()-1].shrink();
    bins.shrink();
    lowerLimits.shrink();
    upperLimits.shrink();


    //
    // Write to cellSets.
    //

    Info<< "Volume bins:" << nl;
    forAll(bins, binI)
    {
        const dynamicLabelList& bin = bins[binI];

        cellSet cells(mesh, "vol" + name(binI), bin.size());

        forAll(bin, i)
        {
            cells.insert(bin[i]);
        }

        Info<< "    " << lowerLimits[binI] << " .. " << upperLimits[binI]
            << "  : writing " << bin.size() << " cells to cellSet "
            << cells.name() << endl;

        cells.write();
    }



    //
    // Convert bins into refinement level.
    //


    // Construct fvMesh to be able to construct volScalarField

    fvMesh fMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        ),
        xferCopy(mesh.points()),   // could we safely re-use the data?
        xferCopy(mesh.faces()),
        xferCopy(mesh.cells())
    );

    // Add the boundary patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> p(patches.size());

    forAll (p, patchI)
    {
        p[patchI] = patches[patchI].clone(fMesh.boundaryMesh()).ptr();
    }

    fMesh.addFvPatches(p);


    // Refinement level
    IOobject refHeader
    (
        "refinementLevel",
        runTime.timeName(),
        polyMesh::defaultRegion,
        runTime
    );

    if (!readLevel && refHeader.headerOk())
    {
        WarningIn(args.executable())
            << "Detected " << refHeader.name() << " file in "
            << polyMesh::defaultRegion <<  " directory. Please remove to"
            << " recreate it or use the -readLevel option to use it"
            << endl;
        return 1;
    }


    labelIOList refLevel
    (
        IOobject
        (
            "refinementLevel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        labelList(mesh.nCells(), 0)
    );

    if (readLevel)
    {
        refLevel = labelIOList(refHeader);
    }

    // Construct volScalarField with same info for post processing
    volScalarField postRefLevel
    (
        IOobject
        (
            "refinementLevel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fMesh,
        dimensionedScalar("zero", dimless/dimTime, 0)
    );

    // Set cell values
    forAll(bins, binI)
    {
        const dynamicLabelList& bin = bins[binI];

        forAll(bin, i)
        {
            refLevel[bin[i]] = bins.size() - binI - 1;
            postRefLevel[bin[i]] = refLevel[bin[i]];
        }
    }


    // For volScalarField: set boundary values to same as cell.
    // Note: could also put
    // zeroGradient b.c. on postRefLevel and do evaluate.
    forAll(postRefLevel.boundaryField(), patchI)
    {
        const polyPatch& pp = patches[patchI];

        fvPatchScalarField& bField = postRefLevel.boundaryField()[patchI];

        Info<< "Setting field for patch "<< endl;

        forAll(bField, faceI)
        {
            label own = mesh.faceOwner()[pp.start() + faceI];

            bField[faceI] = postRefLevel[own];
        }
    }

    Info<< "Determined current refinement level and writing to "
        << postRefLevel.name() << " (as volScalarField; for post processing)"
        << nl
        << polyMesh::defaultRegion/refLevel.name()
        << " (as labelIOList; for meshing)" << nl
        << endl;

    refLevel.write();
    postRefLevel.write();


    // Find out cells to refine to keep to 2:1 refinement level restriction

    // Cells to refine
    cellSet refCells(mesh, "refCells", 100);

    while
    (
        limitRefinementLevel
        (
            mesh,
            refLevel,       // current refinement level
            refCells        // cells to refine
        )
    )
    {}

    if (refCells.size())
    {
        Info<< "Collected " << refCells.size() << " cells that need to be"
            << " refined to get closer to overall 2:1 refinement level limit"
            << nl
            << "Written cells to be refined to cellSet " << refCells.name()
            << nl << endl;

        refCells.write();

        Info<< "After refinement this tool can be run again to see if the 2:1"
            << " limit is observed all over the mesh" << nl << endl;
    }
    else
    {
        Info<< "All cells in the mesh observe the 2:1 refinement level limit"
            << nl << endl;
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
