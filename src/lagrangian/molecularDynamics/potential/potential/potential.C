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

#include "potential.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::potential::setSiteIdList(const IOdictionary& moleculePropertiesDict)
{
    DynamicList<word> siteIdList;

    DynamicList<word> pairPotentialSiteIdList;

    forAll(idList_, i)
    {
        const word& id(idList_[i]);

        if (!moleculePropertiesDict.found(id))
        {
            FatalErrorIn("potential.C") << nl
                << id << " molecule subDict not found"
                << nl << abort(FatalError);
        }

        const dictionary& molDict(moleculePropertiesDict.subDict(id));

        wordList siteIdNames = molDict.lookup("siteIds");

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            if(findIndex(siteIdList, siteId) == -1)
            {
                siteIdList.append(siteId);
            }
        }

        wordList pairPotSiteIds = molDict.lookup("pairPotentialSiteIds");

        forAll(pairPotSiteIds, sI)
        {
            const word& siteId = pairPotSiteIds[sI];

            if(findIndex(siteIdNames, siteId) == -1)
            {
                FatalErrorIn("potential.C") << nl
                    << siteId << " in pairPotentialSiteIds is not in siteIds: "
                    << siteIdNames << nl << abort(FatalError);
            }

            if(findIndex(pairPotentialSiteIdList, siteId) == -1)
            {
                pairPotentialSiteIdList.append(siteId);
            }
        }
    }

    nPairPotIds_ = pairPotentialSiteIdList.size();

    forAll(siteIdList, aSIN)
    {
        const word& siteId = siteIdList[aSIN];

        if(findIndex(pairPotentialSiteIdList, siteId) == -1)
        {
            pairPotentialSiteIdList.append(siteId);
        }
    }

    siteIdList_.transfer(pairPotentialSiteIdList.shrink());
}


void Foam::potential::potential::readPotentialDict()
{
    Info<< nl <<  "Reading potential dictionary:" << endl;

    IOdictionary idListDict
    (
        IOobject
        (
            "idList",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    idList_ = wordList(idListDict.lookup("idList"));

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    setSiteIdList(moleculePropertiesDict);

    wordList pairPotentialSiteIdList
    (
        SubList<word>(siteIdList_, nPairPotIds_)
    );

    Info<< nl << "Unique site ids found: " << siteIdList_
        << nl << "Site Ids requiring a pair potential: "
        << pairPotentialSiteIdList
        << endl;

    wordList tetherSiteIdList(0);

    if (idListDict.found("tetherSiteIdList"))
    {
        tetherSiteIdList = wordList(idListDict.lookup("tetherSiteIdList"));
    }

    IOdictionary potentialDict
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    potentialEnergyLimit_ = readScalar
    (
        potentialDict.lookup("potentialEnergyLimit")
    );

    if (potentialDict.found("removalOrder"))
    {
        wordList remOrd = potentialDict.lookup("removalOrder");

        removalOrder_.setSize(remOrd.size());

        forAll(removalOrder_, rO)
        {
            removalOrder_[rO] = findIndex(idList_, remOrd[rO]);

            if (removalOrder_[rO] == -1)
            {
                FatalErrorIn("potentials.C") << nl
                    << "removalOrder entry: " << remOrd[rO]
                    << " not found in idList."
                    << nl << abort(FatalError);
            }
        }
    }

    // *************************************************************************
    // Pair potentials

    if (!potentialDict.found("pair"))
    {
        FatalErrorIn("potentials.C") << nl
            << "pair potential specification subDict not found"
            << abort(FatalError);
    }

    const dictionary& pairDict = potentialDict.subDict("pair");

    pairPotentials_.buildPotentials
    (
        pairPotentialSiteIdList,
        pairDict,
        mesh_
    );

    // *************************************************************************
    // Tether potentials

    if (tetherSiteIdList.size())
    {
        if (!potentialDict.found("tether"))
        {
            FatalErrorIn("potential.C") << nl
                << "tether potential specification subDict not found"
                << abort(FatalError);
        }

        const dictionary& tetherDict = potentialDict.subDict("tether");

        tetherPotentials_.buildPotentials
        (
            siteIdList_,
            tetherDict,
            tetherSiteIdList
        );
    }

    // *************************************************************************
    // External Forces

    gravity_ = vector::zero;

    if (potentialDict.found("external"))
    {

        Info << nl << "Reading external forces:" << endl;

        const dictionary& externalDict = potentialDict.subDict("external");

        // *********************************************************************
        // gravity

        if (externalDict.found("gravity"))
        {
            gravity_ = externalDict.lookup("gravity");
        }
    }

    Info << nl << tab << "gravity = " << gravity_ << endl;
}


void Foam::potential::potential::readMdInitialiseDict
(
    const IOdictionary& mdInitialiseDict,
    IOdictionary& idListDict
)
{
    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    DynamicList<word> idList;

    DynamicList<word> tetherSiteIdList;

    forAll(mdInitialiseDict.toc(), zone)
    {
        const dictionary& zoneDict = mdInitialiseDict.subDict
        (
            mdInitialiseDict.toc()[zone]
        );

        wordList latticeIds
        (
            zoneDict.lookup("latticeIds")
        );

        forAll(latticeIds, i)
        {
            const word& id = latticeIds[i];

            if(!moleculePropertiesDict.found(id))
            {
                FatalErrorIn("potential.C") << nl
                    << "Molecule type "
                    << id
                    << " not found in moleculeProperties dictionary."
                    << nl
                    << abort(FatalError);
            }

            if (findIndex(idList,id) == -1)
            {
                idList.append(id);
            }
        }

        wordList tetherSiteIds
        (
            zoneDict.lookup("tetherSiteIds")
        );

        forAll(tetherSiteIds, t)
        {
            const word& tetherSiteId = tetherSiteIds[t];

            bool idFound = false;

            forAll(latticeIds, i)
            {
                if (idFound)
                {
                    break;
                }

                const word& id = latticeIds[i];

                wordList siteIds
                (
                    moleculePropertiesDict.subDict(id).lookup("siteIds")
                );

                if (findIndex(siteIds, tetherSiteId) != -1)
                {
                    idFound = true;
                }
            }

            if (idFound)
            {
                tetherSiteIdList.append(tetherSiteId);
            }
            else
            {
                FatalErrorIn("potential.C") << nl
                    << "Tether id  "
                    << tetherSiteId
                    << " not found as a site of any molecule in zone."
                    << nl
                    << abort(FatalError);
            }
        }
    }

    idList_.transfer(idList);

    tetherSiteIdList.shrink();

    idListDict.add("idList", idList_);

    idListDict.add("tetherSiteIdList", tetherSiteIdList);

    setSiteIdList(moleculePropertiesDict);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::potential::potential(const polyMesh& mesh)
:
    mesh_(mesh)
{
    readPotentialDict();
}


Foam::potential::potential
(
    const polyMesh& mesh,
    const IOdictionary& mdInitialiseDict,
    IOdictionary& idListDict
)
:
    mesh_(mesh)
{
    readMdInitialiseDict(mdInitialiseDict, idListDict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::potential::~potential()
{}


// ************************************************************************* //
