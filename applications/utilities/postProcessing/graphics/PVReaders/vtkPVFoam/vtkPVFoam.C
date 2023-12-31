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

#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"

// Foam includes
#include "fvMesh.H"
#include "foamTime.H"
#include "patchZones.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkPVFoam, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "vtkPVFoamAddToSelection.H"
#include "vtkPVFoamUpdateInfoFields.H"

void Foam::vtkPVFoam::resetCounters()
{
    // Reset mesh part ids and sizes
    partInfoVolume_.reset();
    partInfoPatches_.reset();
    partInfoLagrangian_.reset();
    partInfoCellZones_.reset();
    partInfoFaceZones_.reset();
    partInfoPointZones_.reset();
    partInfoCellSets_.reset();
    partInfoFaceSets_.reset();
    partInfoPointSets_.reset();
}


void Foam::vtkPVFoam::reduceMemory()
{
    forAll(regionPolyDecomp_, i)
    {
        regionPolyDecomp_[i].clear();
    }

    forAll(zonePolyDecomp_, i)
    {
        zonePolyDecomp_[i].clear();
    }

    forAll(csetPolyDecomp_, i)
    {
        csetPolyDecomp_[i].clear();
    }

    if (!reader_->GetCacheMesh())
    {
        delete meshPtr_;
        meshPtr_ = nullptr;
    }
}




int Foam::vtkPVFoam::setTime(int nRequest, const double requestTimes[])
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::setTime(";
        for (int requestI = 0; requestI < nRequest; ++requestI)
        {
            if (requestI)
            {
                Info<< ", ";
            }

            Info<< requestTimes[requestI];
        }
        Info << ") - previousIndex = " << timeIndex_ << endl;
    }

    Time& runTime = dbPtr_();

    // Get times list
    instantList Times = runTime.times();

    int nearestIndex = timeIndex_;
    for (int requestI = 0; requestI < nRequest; ++requestI)
    {
        int index = Time::findClosestTimeIndex(Times, requestTimes[requestI]);
        if (index >= 0 && index != timeIndex_)
        {
            nearestIndex = index;
            break;
        }
    }

    if (nearestIndex < 0)
    {
        nearestIndex = 0;
    }


    // see what has changed
    if (timeIndex_ != nearestIndex)
    {
        timeIndex_ = nearestIndex;
        runTime.setTime(Times[nearestIndex], nearestIndex);

        // the fields change each time
        fieldsChanged_ = true;

        if (meshPtr_)
        {
            if (meshPtr_->readUpdate() != polyMesh::UNCHANGED)
            {
                meshChanged_ = true;
            }
        }
        else
        {
            meshChanged_ = true;
        }

        reader_->UpdateProgress(0.05);

        // this seems to be needed for catching Lagrangian fields
        updateInfo();
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::setTime() - selectedTime="
            << Times[nearestIndex].name() << " index=" << timeIndex_
            << "/" << Times.size()
            << " meshChanged=" << Switch(meshChanged_)
            << " fieldsChanged=" << Switch(fieldsChanged_) << endl;
    }

    return nearestIndex;
}


void Foam::vtkPVFoam::updateMeshPartsStatus()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::updateMeshPartsStatus" << endl;
    }

    vtkDataArraySelection* selection = reader_->GetPartSelection();
    label nElem = selection->GetNumberOfArrays();

    if (partStatus_.size() != nElem)
    {
        partStatus_.setSize(nElem);
        partStatus_ = false;
        meshChanged_ = true;
    }

    // this needs fixing if we wish to re-use the datasets
    partDataset_.setSize(nElem);
    partDataset_ = -1;

    // Read the selected mesh parts (zones, patches ...) and add to list
    forAll(partStatus_, partId)
    {
        const int setting = selection->GetArraySetting(partId);

        if (partStatus_[partId] != setting)
        {
            partStatus_[partId] = setting;
            meshChanged_ = true;
        }

        if (debug)
        {
            Info<< "  part[" << partId << "] = "
                << partStatus_[partId]
                << " : " << selection->GetArrayName(partId) << endl;
        }
    }
    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::updateMeshPartsStatus" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPVFoam::vtkPVFoam
(
    const char* const FileName,
    vtkPVFoamReader* reader
)
:
    reader_(reader),
    dbPtr_(nullptr),
    meshPtr_(nullptr),
    meshRegion_(polyMesh::defaultRegion),
    meshDir_(polyMesh::meshSubDir),
    timeIndex_(-1),
    meshChanged_(true),
    fieldsChanged_(true),
    partInfoVolume_("unzoned"),
    partInfoPatches_("patches"),
    partInfoLagrangian_("lagrangian"),
    partInfoCellZones_("cellZone"),
    partInfoFaceZones_("faceZone"),
    partInfoPointZones_("pointZone"),
    partInfoCellSets_("cellSet"),
    partInfoFaceSets_("faceSet"),
    partInfoPointSets_("pointSet")
{
    if (debug)
    {
        Info<< "Foam::vtkPVFoam::vtkPVFoam - " << FileName << endl;
        printMemory();
    }

    // avoid argList and get rootPath/caseName directly from the file
    fileName fullCasePath(fileName(FileName).path());

    if (!isDir(fullCasePath))
    {
        return;
    }
    if (fullCasePath == ".")
    {
        fullCasePath = cwd();
    }

    // Set the case as an environment variable - some BCs might use this
    if (fullCasePath.name().find("processor", 0) == 0)
    {
        const fileName globalCase = fullCasePath.path();

        setEnv("FOAM_CASE", globalCase, true);
        setEnv("FOAM_CASENAME", globalCase.name(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
        setEnv("FOAM_CASENAME", fullCasePath.name(), true);
    }

    // look for 'case{region}.FOAM'
    // could be stringent and insist the prefix match the directory name...
    // Note: cannot use fileName::name() due to the embedded '{}'
    string caseName(fileName(FileName).lessExt());
    string::size_type beg = caseName.find_last_of("/{");
    string::size_type end = caseName.find('}', beg);

    if
    (
        beg != string::npos && caseName[beg] == '{'
     && end != string::npos && end == caseName.size()-1
    )
    {
        meshRegion_ = caseName.substr(beg+1, end-beg-1);

        // some safety
        if (!meshRegion_.size())
        {
            meshRegion_ = polyMesh::defaultRegion;
        }

        if (meshRegion_ != polyMesh::defaultRegion)
        {
            meshDir_ = meshRegion_/polyMesh::meshSubDir;
        }
    }

    if (debug)
    {
        Info<< "fullCasePath=" << fullCasePath << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
            << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
            << "region=" << meshRegion_ << endl;
    }

    // Create time object
    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            fileName(fullCasePath.path()),
            fileName(fullCasePath.name())
        )
    );

    dbPtr_().functionObjects().off();

    updateInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPVFoam::~vtkPVFoam()
{
    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::~vtkPVFoam" << endl;
    }

    delete meshPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::updateInfo()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::updateInfo"
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "] timeIndex="
            << timeIndex_ << endl;
    }

    resetCounters();

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();

    // enable 'internalMesh' on the first call
    // or preserve the enabled selections
    stringList enabledEntries;
    if (!partSelection->GetNumberOfArrays() && !meshPtr_)
    {
        enabledEntries.setSize(1);
        enabledEntries[0] = "internalMesh";
    }
    else
    {
        enabledEntries = getSelectedArrayEntries(partSelection);
    }

    // Clear current mesh parts list
    partSelection->RemoveAllArrays();

    // Update mesh parts list - add Lagrangian at the bottom
    updateInfoInternalMesh();
    updateInfoPatches();
    updateInfoSets();
    updateInfoZones();
    updateInfoLagrangian();

    // restore the enabled selections
    setSelectedArrayEntries(partSelection, enabledEntries);

    if (meshChanged_)
    {
        fieldsChanged_ = true;
    }

    // Update volume, point and lagrangian fields
    updateInfoFields<fvPatchField, volMesh>
    (
        reader_->GetVolFieldSelection()
    );
    updateInfoFields<pointPatchField, pointMesh>
    (
        reader_->GetPointFieldSelection()
    );
    updateInfoLagrangianFields();

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::updateInfo" << endl;
    }

}


void Foam::vtkPVFoam::updateFoamMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::updateFoamMesh" << endl;
        printMemory();
    }

    if (!reader_->GetCacheMesh())
    {
        delete meshPtr_;
        meshPtr_ = nullptr;
    }

    // Check to see if the FOAM mesh has been created
    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Creating Foam mesh for region " << meshRegion_
                << " at time=" << dbPtr_().timeName()
                << endl;

        }

        meshPtr_ = new fvMesh
        (
            IOobject
            (
                meshRegion_,
                dbPtr_().timeName(),
                dbPtr_(),
                IOobject::MUST_READ
            )
        );

        meshChanged_ = true;
    }
    else
    {
        if (debug)
        {
            Info<< "Using existing Foam mesh" << endl;
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::updateFoamMesh" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::Update
(
    vtkMultiBlockDataSet* output,
    vtkMultiBlockDataSet* lagrangianOutput
)
{
    if (debug)
    {
        cout<< "<beg> Foam::vtkPVFoam::Update - output with "
            << output->GetNumberOfBlocks() << " and "
            << lagrangianOutput->GetNumberOfBlocks() << " blocks\n";
        output->Print(cout);
        lagrangianOutput->Print(cout);
        printMemory();
    }
    reader_->UpdateProgress(0.1);

    // Set up mesh parts selection(s)
    updateMeshPartsStatus();

    reader_->UpdateProgress(0.15);

    // Update the Foam mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.4);

    // Convert meshes - start port0 at block=0
    int blockNo = 0;

    convertMeshVolume(output, blockNo);
    convertMeshPatches(output, blockNo);
    reader_->UpdateProgress(0.6);

    if (reader_->GetIncludeZones())
    {
        convertMeshCellZones(output, blockNo);
        convertMeshFaceZones(output, blockNo);
        convertMeshPointZones(output, blockNo);
        reader_->UpdateProgress(0.65);
    }

    if (reader_->GetIncludeSets())
    {
        convertMeshCellSets(output, blockNo);
        convertMeshFaceSets(output, blockNo);
        convertMeshPointSets(output, blockNo);
        reader_->UpdateProgress(0.7);
    }

#ifdef VTKPVFOAM_DUALPORT
    // restart port1 at block=0
    blockNo = 0;
#endif
    convertMeshLagrangian(lagrangianOutput, blockNo);

    reader_->UpdateProgress(0.8);

    // Update fields
    convertVolFields(output);
    convertPointFields(output);
    convertLagrangianFields(lagrangianOutput);
    reader_->UpdateProgress(0.95);

    meshChanged_ = fieldsChanged_ = false;
}


void Foam::vtkPVFoam::CleanUp()
{
    // reclaim some memory
    reduceMemory();
    reader_->UpdateProgress(1.0);
}


double* Foam::vtkPVFoam::findTimes(int& nTimeSteps)
{
    int nTimes = 0;
    double* tsteps = nullptr;

    if (dbPtr_.valid())
    {
        Time& runTime = dbPtr_();
        instantList timeLst = runTime.times();

        // find the first time for which this mesh appears to exist
        label timeI = 0;
        for (; timeI < timeLst.size(); ++timeI)
        {
            const word& timeName = timeLst[timeI].name();

            if
            (
                isFile(runTime.path()/timeName/meshDir_/"points")
             && IOobject("points", timeName, meshDir_, runTime).headerOk()
            )
            {
                break;
            }
        }

        nTimes = timeLst.size() - timeI;

        // always skip "constant" time if possible
        if (timeI == 0 && nTimes > 1)
        {
            timeI = 1;
            --nTimes;
        }

        if (nTimes)
        {
            tsteps = new double[nTimes];
            for (label stepI = 0; stepI < nTimes; ++stepI, ++timeI)
            {
                tsteps[stepI] = timeLst[timeI].value();
            }
        }
    }
    else
    {
        if (debug)
        {
            cout<< "no valid dbPtr:\n";
        }
    }

    // vector length returned via the parameter
    nTimeSteps = nTimes;

    return tsteps;
}


void Foam::vtkPVFoam::addPatchNames(vtkRenderer* renderer)
{
    // Remove any patch names previously added to the renderer
    removePatchNames(renderer);

    // get the display patches, strip off any suffix
    wordHashSet selectedPatches = getSelected
    (
        reader_->GetPartSelection(),
        partInfoPatches_
    );

    if (!selectedPatches.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::addPatchNames" << nl
            <<"... add patches: " << selectedPatches << endl;
    }

    const polyBoundaryMesh& pbMesh = meshPtr_->boundaryMesh();

    // Find the total number of zones
    // Each zone will take the patch name
    // Number of zones per patch ... zero zones should be skipped
    labelList nZones(pbMesh.size(), 0);

    // Per global zone number the average face centre position
    DynamicList<point> zoneCentre(pbMesh.size());

    if (debug)
    {
        Info<< "... determining patch zones" << endl;
    }

    // Loop through all patches to determine zones, and centre of each zone
    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only include the patch if it is selected
        if (!selectedPatches.found(pp.name()))
        {
            continue;
        }

        const labelListList& edgeFaces = pp.edgeFaces();
        const vectorField& n = pp.faceNormals();

        boolList featEdge(pp.nEdges(), false);

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = edgeFaces[edgeI];

            if (eFaces.size() == 1)
            {
                // Note: could also do ones with > 2 faces but this gives
                // too many zones for baffles
                featEdge[edgeI] = true;
            }
            else if (mag(n[eFaces[0]] & n[eFaces[1]]) < 0.5)
            {
                featEdge[edgeI] = true;
            }
        }

        // Do topological analysis of patch, find disconnected regions
        patchZones pZones(pp, featEdge);

        nZones[patchI] = pZones.nZones();

        labelList zoneNFaces(pZones.nZones(), 0);

        // Save start of information for current patch
        label patchStart = zoneCentre.size();

        // Create storage for additional zone centres
        forAll(zoneNFaces, zoneI)
        {
            zoneCentre.append(vector::zero);
        }

        // Do averaging per individual zone
        forAll(pp, faceI)
        {
            label zoneI = pZones[faceI];
            zoneCentre[patchStart+zoneI] += pp[faceI].centre(pp.points());
            zoneNFaces[zoneI]++;
        }

        for (label i=0; i<nZones[patchI]; i++)
        {
            zoneCentre[patchStart + i] /= zoneNFaces[i];
        }
    }

    // Count number of zones we're actually going to display. This is truncated
    // to a max per patch

    const label MAXPATCHZONES = 20;

    label displayZoneI = 0;

    forAll(pbMesh, patchI)
    {
        displayZoneI += min(MAXPATCHZONES, nZones[patchI]);
    }


    zoneCentre.shrink();

    if (debug)
    {
        Info<< "patch zone centres = " << zoneCentre << nl
            << "displayed zone centres = " << displayZoneI << nl
            << "zones per patch = " << nZones << endl;
    }

    // Set the size of the patch labels to max number of zones
    patchTextActorsPtrs_.setSize(displayZoneI);

    if (debug)
    {
        Info<< "constructing patch labels" << endl;
    }

    // Actor index
    displayZoneI = 0;

    // Index in zone centres
    label globalZoneI = 0;

    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only selected patches will have a non-zero number of zones
        label nDisplayZones = min(MAXPATCHZONES, nZones[patchI]);
        label increment = 1;
        if (nZones[patchI] >= MAXPATCHZONES)
        {
            increment = nZones[patchI]/MAXPATCHZONES;
        }

        for (label i = 0; i < nDisplayZones; i++)
        {
            if (debug)
            {
                Info<< "patch name = " << pp.name() << nl
                    << "anchor = " << zoneCentre[globalZoneI] << nl
                    << "globalZoneI = " << globalZoneI << endl;
            }

            vtkTextActor* txt = vtkTextActor::New();

            txt->SetInput(pp.name().c_str());

            // Set text properties
            vtkTextProperty* tprop = txt->GetTextProperty();
            tprop->SetFontFamilyToArial();
            tprop->BoldOff();
            tprop->ShadowOff();
            tprop->SetLineSpacing(1.0);
            tprop->SetFontSize(12);
            tprop->SetColor(1.0, 0.0, 0.0);
            tprop->SetJustificationToCentered();

            // Set text to use 3-D world co-ordinates
            txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();

            txt->GetPositionCoordinate()->SetValue
            (
                zoneCentre[globalZoneI].x(),
                zoneCentre[globalZoneI].y(),
                zoneCentre[globalZoneI].z()
            );

            // Add text to each renderer
            renderer->AddViewProp(txt);

            // Maintain a list of text labels added so that they can be
            // removed later
            patchTextActorsPtrs_[displayZoneI] = txt;

            globalZoneI += increment;
            displayZoneI++;
        }
    }

    // Resize the patch names list to the actual number of patch names added
    patchTextActorsPtrs_.setSize(displayZoneI);

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::addPatchNames" << endl;
    }
}


void Foam::vtkPVFoam::removePatchNames(vtkRenderer* renderer)
{
    forAll(patchTextActorsPtrs_, patchI)
    {
        renderer->RemoveViewProp(patchTextActorsPtrs_[patchI]);
        patchTextActorsPtrs_[patchI]->Delete();
    }
    patchTextActorsPtrs_.clear();
}


void Foam::vtkPVFoam::PrintSelf(ostream& os, vtkIndent indent) const
{
    os  << indent << "Number of nodes: "
        << (meshPtr_ ? meshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (meshPtr_ ? meshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_.valid() ? dbPtr_().times().size() : 0) << endl;

    os  << indent << "mesh region: " << meshRegion_ << "\n";
}

// ************************************************************************* //
