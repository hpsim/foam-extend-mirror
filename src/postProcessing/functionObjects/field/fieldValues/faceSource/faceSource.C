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

#include "faceSource.H"
#include "fvMesh.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "sampledSurface.H"
#include "mergePoints.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<fieldValues::faceSource::sourceType, 3>::names[] =
    {
        "faceZone",
        "patch",
        "sampledSurface"
    };


    template<>
    const char* NamedEnum<fieldValues::faceSource::operationType, 17>::names[] =
    {
        "none",
        "sum",
        "weightedSum",
        "sumMag",
        "weightedSumMag",
        "sumDirection",
        "sumDirectionBalance",
        "average",
        "weightedAverage",
        "areaAverage",
        "weightedAreaAverage",
        "areaIntegrate",
        "min",
        "max",
        "CoV",
        "areaNormalAverage",
        "areaNormalIntegrate"
    };

    namespace fieldValues
    {
        defineTypeNameAndDebug(faceSource, 0);
        addToRunTimeSelectionTable(fieldValue, faceSource, dictionary);
    }
}


const Foam::NamedEnum<Foam::fieldValues::faceSource::sourceType, 3>
    Foam::fieldValues::faceSource::sourceTypeNames_;

const Foam::NamedEnum<Foam::fieldValues::faceSource::operationType, 17>
    Foam::fieldValues::faceSource::operationTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSource::setFaceZoneFaces()
{
    label zoneId = mesh().faceZones().findZoneID(sourceName_);

    if (zoneId < 0)
    {
        FatalErrorIn("faceSource::faceSource::setFaceZoneFaces()")
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Unknown face zone name: " << sourceName_
            << ". Valid face zones are: " << mesh().faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh().faceZones()[zoneId];

    dynamicLabelList faceIds(fZone.size());
    dynamicLabelList facePatchIds(fZone.size());
    dynamicLabelList faceSigns(fZone.size());

    forAll(fZone, i)
    {
        label faceI = fZone[i];

        label faceId = -1;
        label facePatchId = -1;
        if (mesh().isInternalFace(faceI))
        {
            faceId = faceI;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh().boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh().boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                //HR: need to review!
                //if (refCast<const coupledPolyPatch>(pp).owner())
                //{
                    faceId = pp.whichFace(faceI);
                //}
                //else
                //{
                //faceId = -1;
                //}
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = faceI - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSigns.append(-1);
            }
            else
            {
                faceSigns.append(1);
            }
            faceIds.append(faceId);
            facePatchIds.append(facePatchId);
        }
    }

    faceId_.transfer(faceIds);
    facePatchId_.transfer(facePatchIds);
    faceSign_.transfer(faceSigns);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    if (debug)
    {
        Pout<< "Original face zone size = " << fZone.size()
            << ", new size = " << faceId_.size() << endl;
    }
}


void Foam::fieldValues::faceSource::setPatchFaces()
{
    const label patchId = mesh().boundaryMesh().findPatchID(sourceName_);

    if (patchId < 0)
    {
        FatalErrorIn("faceSource::constructFaceAddressing()")
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Unknown patch name: " << sourceName_
            << ". Valid patch names are: "
            << mesh().boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const polyPatch& pp = mesh().boundaryMesh()[patchId];

    label nFaces = pp.size();
    if (isA<emptyPolyPatch>(pp))
    {
        nFaces = 0;
    }

    faceId_.setSize(nFaces);
    facePatchId_.setSize(nFaces);
    faceSign_.setSize(nFaces);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    forAll(faceId_, faceI)
    {
        faceId_[faceI] = faceI;
        facePatchId_[faceI] = patchId;
        faceSign_[faceI] = 1;
    }
}


void Foam::fieldValues::faceSource::sampledSurfaceFaces(const dictionary& dict)
{
    surfacePtr_ = sampledSurface::New
    (
        name_,
        mesh(),
        dict.subDict("sampledSurfaceDict")
    );
    surfacePtr_().update();
    nFaces_ = returnReduce(surfacePtr_().faces().size(), sumOp<label>());
}


void Foam::fieldValues::faceSource::combineMeshGeometry
(
    faceList& faces,
    pointField& points
) const
{
    List<faceList> allFaces(Pstream::nProcs());
    List<pointField> allPoints(Pstream::nProcs());

    labelList globalFacesIs(faceId_);
    forAll(globalFacesIs, i)
    {
        if (facePatchId_[i] != -1)
        {
            label patchI = facePatchId_[i];
            globalFacesIs[i] += mesh().boundaryMesh()[patchI].start();
        }
    }

    // Add local faces and points to the all* lists
    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh().faces(), globalFacesIs),
        mesh().points()
    );
    allFaces[Pstream::myProcNo()] = pp.localFaces();
    allPoints[Pstream::myProcNo()] = pp.localPoints();

    Pstream::gatherList(allFaces);
    Pstream::gatherList(allPoints);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, procI)
    {
        nFaces += allFaces[procI].size();
        nPoints += allPoints[procI].size();
    }

    faces.setSize(nFaces);
    points.setSize(nPoints);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const faceList& fcs = allFaces[Pstream::myProcNo()];
        forAll(fcs, i)
        {
            const face& f = fcs[i];
            face& newF = faces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        forAll(pts, i)
        {
            points[nPoints++] = pts[i];
        }
    }

    // Other proc data follows
    forAll(allFaces, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            const faceList& fcs = allFaces[procI];
            forAll(fcs, i)
            {
                const face& f = fcs[i];
                face& newF = faces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[procI];
            forAll(pts, i)
            {
                points[nPoints++] = pts[i];
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newPoints;
    bool hasMerged = mergePoints
    (
        points,
        SMALL,
        false,
        oldToNew,
        newPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << points.size()
                << " down to " << newPoints.size() << " points" << endl;
        }

        points.transfer(newPoints);
        forAll(faces, i)
        {
            inplaceRenumber(oldToNew, faces[i]);
        }
    }
}


void Foam::fieldValues::faceSource::combineSurfaceGeometry
(
    faceList& faces,
    pointField& points
) const
{
    if (surfacePtr_.valid())
    {
        const sampledSurface& s = surfacePtr_();

        if (Pstream::parRun())
        {
            // dimension as fraction of mesh bounding box
            scalar mergeDim = 1e-10*mesh().bounds().mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
}


Foam::scalar Foam::fieldValues::faceSource::totalArea() const
{
    scalar totalArea;

    if (surfacePtr_.valid())
    {
        totalArea = gSum(surfacePtr_().magSf());
    }
    else
    {
        totalArea = gSum(filterField(mesh().magSf(), false));
    }

    return totalArea;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSource::initialise(const dictionary& dict)
{
    dict.lookup("sourceName") >> sourceName_;

    switch (source_)
    {
        case stFaceZone:
        {
            setFaceZoneFaces();
            break;
        }
        case stPatch:
        {
            setPatchFaces();
            break;
        }
        case stSampledSurface:
        {
            sampledSurfaceFaces(dict);
            break;
        }
        default:
        {
            FatalErrorIn("faceSource::initialise()")
                << type() << " " << name_ << ": "
                << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                << nl << "    Unknown source type. Valid source types are:"
                << sourceTypeNames_.sortedToc() << nl << exit(FatalError);
        }
    }

    if (nFaces_ == 0)
    {
        WarningIn
        (
            "Foam::fieldValues::faceSource::initialise(const dictionary&)"
        )
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Source has no faces - deactivating" << endl;

        active_ = false;
        return;
    }

    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    totalArea_ = totalArea();

    Info<< type() << " " << name_ << ":" << nl
        << "    total faces  = " << nFaces_
        << nl
        << "    total area   = " << totalArea_
        << nl;

    if (dict.readIfPresent("weightField", weightFieldName_))
    {
        Info<< "    weight field = " << weightFieldName_ << nl;

        if (source_ == stSampledSurface)
        {
            FatalIOErrorIn
            (
                "void Foam::fieldValues::faceSource::initialise"
                "("
                    "const dictionary&"
                ")",
                dict
            )
                << "Cannot use weightField for a sampledSurface"
                << exit(FatalIOError);
        }
    }

    if (dict.found("orientedWeightField"))
    {
        if (weightFieldName_ == "none")
        {
            dict.lookup("orientedWeightField") >>  weightFieldName_;
            Info<< "    weight field = " << weightFieldName_ << nl;
            orientWeightField_ = true;
        }
        else
        {
            FatalIOErrorIn
            (
                "void Foam::fieldValues::faceSource::initialise"
                "("
                    "const dictionary&"
                ")",
                dict
            )
                << "Either weightField or orientedWeightField can be supplied, "
                << "but not both"
                << exit(FatalIOError);
        }
    }

    wordList orientedFields;
    if (dict.readIfPresent("orientedFields", orientedFields))
    {
        orientedFieldsStart_ = fields_.size();
        fields_.append(orientedFields);
    }

    if (dict.readIfPresent("scaleFactor", scaleFactor_))
    {
        Info<< "    scale factor = " << scaleFactor_ << nl;
    }

    Info<< nl << endl;

    if (valueOutput_)
    {
        const word surfaceFormat(dict.lookup("surfaceFormat"));

        surfaceWriterPtr_.reset
        (
            surfaceWriter::New
            (
                surfaceFormat,
                dict.subOrEmptyDict("formatOptions").
                    subOrEmptyDict(surfaceFormat)
            ).ptr()
        );
    }
}


void Foam::fieldValues::faceSource::writeFileHeader(const label i)
{
    writeCommented(file(), "Source : ");
    file() << sourceTypeNames_[source_] << " " << sourceName_ << endl;
    writeCommented(file(), "Faces  : ");
    file() << nFaces_ << endl;
    writeCommented(file(), "Area   : ");
    file() << totalArea_ << endl;

    writeCommented(file(), "Time");
    if (writeArea_)
    {
        file() << tab << "Area";
    }

    forAll(fields_, i)
    {
        file()
            << tab << operationTypeNames_[operation_]
            << "(" << fields_[i] << ")";
    }

    file() << endl;
}


template<>
Foam::scalar Foam::fieldValues::faceSource::processValues
(
    const Field<scalar>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            vector n(dict_.lookup("direction"));
            return sum(pos(values*(Sf & n))*mag(values));
        }
        case opSumDirectionBalance:
        {
            vector n(dict_.lookup("direction"));
            const scalarField nv(values*(Sf & n));

            return sum(pos(nv)*mag(values) - neg(nv)*mag(values));
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


template<>
Foam::vector Foam::fieldValues::faceSource::processValues
(
    const Field<vector>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            vector n(dict_.lookup("direction"));
            n /= mag(n) + ROOTVSMALL;
            const scalarField nv(n & values);

            return sum(pos(nv)*n*(nv));
        }
        case opSumDirectionBalance:
        {
            vector n(dict_.lookup("direction"));
            n /= mag(n) + ROOTVSMALL;
            const scalarField nv(n & values);

            return sum(pos(nv)*n*(nv));
        }
        case opAreaNormalAverage:
        {
            scalar result = sum(values & Sf)/sum(mag(Sf));
            return vector(result, 0.0, 0.0);
        }
        case opAreaNormalIntegrate:
        {
            scalar result = sum(values & Sf);
            return vector(result, 0.0, 0.0);
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValues::faceSource::faceSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    fieldValue(name, obr, dict, typeName, loadFromFiles),
    surfaceWriterPtr_(nullptr),
    source_(sourceTypeNames_.read(dict.lookup("source"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldName_("none"),
    orientWeightField_(false),
    orientedFieldsStart_(labelMax),
    scaleFactor_(1.0),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::faceSource::~faceSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::faceSource::read(const dictionary& dict)
{
    fieldValue::read(dict);

    if (active_)
    {
        initialise(dict);
    }
}


void Foam::fieldValues::faceSource::write()
{
    fieldValue::write();

    if (active_)
    {
        if (surfacePtr_.valid())
        {
            surfacePtr_().update();
        }

        if (Pstream::master())
        {
            file() << obr_.time().value();
        }

        if (writeArea_)
        {
            totalArea_ = totalArea();
            if (Pstream::master())
            {
                file() << tab << totalArea_;
            }
            if (log_) Info<< "    total area = " << totalArea_ << endl;
        }

        // construct weight field. Note: zero size means weight = 1
        scalarField weightField;
        if (weightFieldName_ != "none")
        {
            weightField =
                getFieldValues<scalar>
                (
                    weightFieldName_,
                    true,
                    orientWeightField_
                );
        }

        // Combine onto master
        combineFields(weightField);

        // process the fields
        forAll(fields_, i)
        {
            const word& fieldName = fields_[i];
            bool ok = false;

            bool orient = i >= orientedFieldsStart_;
            ok = ok || writeValues<scalar>(fieldName, weightField, orient);
            ok = ok || writeValues<vector>(fieldName, weightField, orient);
            ok = ok
              || writeValues<sphericalTensor>(fieldName, weightField, orient);
            ok = ok || writeValues<symmTensor>(fieldName, weightField, orient);
            ok = ok || writeValues<tensor>(fieldName, weightField, orient);

            if (!ok)
            {
                WarningIn("void Foam::fieldValues::faceSource::write()")
                    << "Requested field " << fieldName
                    << " not found in database and not processed"
                    << endl;
            }
        }

        if (Pstream::master())
        {
            file()<< endl;
        }

        if (log_) Info<< endl;
    }
}


// ************************************************************************* //
