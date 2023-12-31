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

#include "linearValveFvMesh.H"
#include "foamTime.H"
#include "slidingInterface.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearValveFvMesh, 0);

    addToRunTimeSelectionTable(topoChangerFvMesh, linearValveFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearValveFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void linearValveFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void linearValveFvMesh::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(1);

    // Add an empty zone for cut points

    pz[0] = new pointZone
    (
        "cutPointZone",
        labelList(0),
        0,
        pointZones()
    );


    // Do face zones for slider

    List<faceZone*> fz(3);

    // Inner slider
    const word innerSliderName(motionDict_.subDict("slider").lookup("inside"));
    const polyPatch& innerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(innerSliderName)];

    labelList isf(innerSlider.size());

    forAll (isf, i)
    {
        isf[i] = innerSlider.start() + i;
    }

    fz[0] = new faceZone
    (
        "insideSliderZone",
        isf,
        boolList(innerSlider.size(), false),
        0,
        faceZones()
    );

    // Outer slider
    const word outerSliderName(motionDict_.subDict("slider").lookup("outside"));
    const polyPatch& outerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(outerSliderName)];

    labelList osf(outerSlider.size());

    forAll (osf, i)
    {
        osf[i] = outerSlider.start() + i;
    }

    fz[1] = new faceZone
    (
        "outsideSliderZone",
        osf,
        boolList(outerSlider.size(), false),
        1,
        faceZones()
    );

    // Add empty zone for cut faces
    fz[2] = new faceZone
    (
        "cutFaceZone",
        labelList(0),
        boolList(0, false),
        2,
        faceZones()
    );

    List<cellZone*> cz(0);

    Info << "Adding point, face and cell zones" << endl;
    addZones(pz, fz, cz);

    // Add a topology modifier
    Info << "Adding topology modifiers" << endl;
    topoChanger_.setSize(1);
    topoChanger_.set
    (
        0,
        new slidingInterface
        (
            "mixerSlider",
            0,
            topoChanger_,
            outerSliderName + "Zone",
            innerSliderName + "Zone",
            "cutPointZone",
            "cutFaceZone",
            outerSliderName,
            innerSliderName,
            slidingInterface::INTEGRAL,   // Edge matching algorithm
            true,                         // Attach-detach action
            intersection::VISIBLE         // Projection algorithm
        )
    );

    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


void Foam::linearValveFvMesh::deactivateSliders()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else
        {
            FatalErrorIn("void Foam::linearValveFvMesh::deactivateSliders()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::linearValveFvMesh::activateSliders()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable sliding interface
    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorIn("void Foam::linearValveFvMesh::activateSliders()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


bool Foam::linearValveFvMesh::attached() const
{
    const polyTopoChanger& topoChanges = topoChanger_;

    bool result = false;

    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanges[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            if
            (
                result
             != refCast<const slidingInterface>(topoChanges[modI]).attached()
            )
            {
                FatalErrorIn("bool linearValveFvMesh::attached() const")
                    << "Slider " << modI << " named "
                    << topoChanges[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    if (result)
    {
        Info << "linearValveFvMesh: attached!" << endl;
    }
    else
    {
        Info << "linearValveFvMesh: detached!" << endl;
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearValveFvMesh::linearValveFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    msPtr_(motionSolver::New(*this))
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearValveFvMesh::~linearValveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::linearValveFvMesh::update()
{
    // Detaching the interface
    if (attached())
    {
        Info << "Decoupling sliding interfaces" << endl;
        activateSliders();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

        if (topoChangeMap1->morphing())
        {
            msPtr_->updateMesh(topoChangeMap1());
        }
    }
    else
    {
        Info << "Sliding interfaces decoupled" << endl;
    }

    // Perform mesh motion
    deactivateSliders();

    // Changing topology by hand
    {
        autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

        if (topoChangeMap2->morphing())
        {
            msPtr_->updateMesh(topoChangeMap2());

            if (topoChangeMap2->hasMotionPoints())
            {
                Info << "Topology change; executing pre-motion" << endl;
                movePoints(topoChangeMap2->preMotionPoints());
            }
        }
    }

    // Solve for motion
    msPtr_->solve();

    movePoints(msPtr_->curPoints());

    // Attach the interface
    Info << "Coupling sliding interfaces" << endl;
    activateSliders();

    // Changing topology by hand
    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        Info << "Moving points post slider attach" << endl;

        bool localMorphing3 = topoChangeMap3->morphing();
        bool globalMorphing3 = localMorphing3;

        reduce(globalMorphing3, orOp<bool>());

        if (globalMorphing3)
        {
            pointField newPoints = allPoints();

            if (localMorphing3)
            {
                msPtr_->updateMesh(topoChangeMap3());

                pointField mappedOldPointsNew(newPoints.size());

                mappedOldPointsNew.map
                (
                    oldPointsNew,
                    topoChangeMap3->pointMap()
                );

                // Solve the correct mesh motion to make sure motion fluxes
                // are solved for and not mapped
                // Note: using setOldPoints instead of movePoints.
                // HJ, 23/Aug/2015
                setOldPoints(mappedOldPointsNew);

                resetMotion();
                setV0();

                fvMesh::movePoints(newPoints);
            }
            else
            {
                // No local topological change.  Execute double motion for
                // sync with topological changes
                // Note: using setOldPoints instead of movePoints.
                // HJ, 23/Aug/2015
                setOldPoints(oldPointsNew);

                resetMotion();
                setV0();

                // Set new point motion
                fvMesh::movePoints(newPoints);
            }
        }
    }

    Info << "Sliding interfaces coupled: " << attached() << endl;

    return true;
}


// ************************************************************************* //
