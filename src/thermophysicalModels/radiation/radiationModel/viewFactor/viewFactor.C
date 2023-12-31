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

#include "viewFactor.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"
#include "greyDiffusiveViewFactorFixedValueFvPatchScalarField.H"
#include "typeInfo.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(viewFactor, 0);
        addToRadiationRunTimeSelectionTables(viewFactor);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::viewFactor::initialise()
{
    const polyBoundaryMesh& coarsePatches = coarseMesh_.boundaryMesh();
    const volScalarField::GeometricBoundaryField& Qrp = Qr_.boundaryField();

    label count = 0;
    forAll (Qrp, patchI)
    {
        const fvPatchScalarField& QrPatchI = Qrp[patchI];

        if ((isA<fixedValueFvPatchScalarField>(QrPatchI)))
        {
            selectedPatches_[count] = QrPatchI.patch().index();
            nLocalCoarseFaces_ += coarsePatches[patchI].size();
            count++;
        }
    }

    selectedPatches_.resize(count);

    if (debug)
    {
        Pout<< "Selected patches:" << selectedPatches_ << endl;
        Pout<< "Number of coarse faces:" << nLocalCoarseFaces_ << endl;
    }

    totalNCoarseFaces_ = nLocalCoarseFaces_;
    reduce(totalNCoarseFaces_, sumOp<label>());

    if (Pstream::master())
    {
        Info<< "Total number of clusters : " << totalNCoarseFaces_ << endl;
    }

    labelListIOList subMap
    (
        IOobject
        (
            "subMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList constructMap
    (
        IOobject
        (
            "constructMap",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelIOList consMapDim
    (
        IOobject
        (
            "constructMapDim",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (consMapDim.empty())
    {
        FatalErrorIn("void radiation::viewFactor::initialise()")
            << "constructMap/consMapDim is empty: radiation model cannot be "
            << "constructed"
            << abort(FatalError);
    }

    map_.reset
    (
        new mapDistribute
        (
            consMapDim[0],
            subMap,
            constructMap
        )
    );

    scalarListIOList FmyProc
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelListIOList globalFaceFaces
    (
        IOobject
        (
            "globalFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    List<labelListList> globalFaceFacesProc(Pstream::nProcs());
    globalFaceFacesProc[Pstream::myProcNo()] = globalFaceFaces;
    Pstream::gatherList(globalFaceFacesProc);

    List<scalarListList> F(Pstream::nProcs());
    F[Pstream::myProcNo()] = FmyProc;
    Pstream::gatherList(F);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Calculate global face areas
    scalarField compactCoarseMagSf(map_->constructSize(), 0.0);
    DynamicList<scalar> localCoarseMagSf(nLocalCoarseFaces_);

    forAll (selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID];
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID];

        scalarList magSf(pp.size(), 0.0);

        if (pp.size() > 0)
        {
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom) + 1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            forAll (coarseToFine, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                UIndirectList<scalar> fineSf
                (
                    sf,
                    fineFaces
                );
                magSf[coarseI] = sum(fineSf());
            }
        }

        localCoarseMagSf.append(magSf);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseMagSf, nLocalCoarseFaces_).assign(localCoarseMagSf);

    // Distribute data
    map_->distribute(compactCoarseMagSf);

    // Distribute local global ID
    labelList compactGlobalIds(map_->constructSize(), 0.0);

    labelList localGlobalIds(nLocalCoarseFaces_);

    for(label k = 0; k < nLocalCoarseFaces_; k++)
    {
        localGlobalIds[k] = globalNumbering.toGlobal(Pstream::myProcNo(), k);
    }

    SubList<label>
    (
        compactGlobalIds,
        nLocalCoarseFaces_
    ).assign(localGlobalIds);

    map_->distribute(compactGlobalIds);

    // Create global size vectors
    scalarField coarseMagSf(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll (compactCoarseMagSf, i)
    {
        coarseMagSf[compactGlobalIds[i]] = compactCoarseMagSf[i];
    }

    Pstream::listCombineGather(coarseMagSf, maxEqOp<scalar>());

    Pstream::listCombineScatter(coarseMagSf);


    if (Pstream::master())
    {
        Fmatrix_.reset
        (
            new scalarSquareMatrix(totalNCoarseFaces_, 0.0)
        );

        Info<< "Insert elements in the matrix..." << endl;

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            insertMatrixElements
            (
                globalNumbering,
                procI,
                globalFaceFacesProc[procI],
                F[procI],
                Fmatrix_()
            );
        }


        bool smoothing = readBool(coeffs_.lookup("smoothing"));
        if (smoothing)
        {
            Info<< "Smoothing the matrix..." << endl;

            // HR 19.02.19: Old smoothing messed up the reciprocity condition:
            // A_i*F_ij = A_j*F_ji. The new algorithm distributes the error
            // in two passes moving downwards and upwards through the rows.
            // On the way down, only the upper triangular is updated from the
            // closedness condition: sum_j(F_ij) = 1; the lower triangular is
            // updated from reciprocity. The reverse happen on the way up.

            if(debug)
            {
                scalar maxErr = -GREAT;
                scalar minSumF = GREAT;
                scalar maxSumF = -GREAT;
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    scalar sumF = 0.0;
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    {
                        sumF += Fmatrix_()[i][j];
                    }
                    maxErr = max(maxErr, mag(sumF - 1.0));
                    minSumF = min(minSumF, sumF);
                    maxSumF = max(maxSumF, sumF);
                }

                Info << "max err = " << maxErr
                    << " min sumF = " << minSumF
                    << " max sumF = " << maxSumF << endl;
            }

            label iter = 0;
            label maxIter = 20;
            while(true)
            {
                iter++;
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    scalar sumF1 = 0.0;
                    for (label j=0; j<i; j++)
                    {
                        sumF1 += Fmatrix_()[i][j];
                    }

                    if(sumF1 > 1.0)
                    {
                        continue;
                    }

                    scalar sumF2 = 0.0;
                    for (label j=i+1; j<totalNCoarseFaces_; j++)
                    {
                        sumF2 += Fmatrix_()[i][j];
                    }

                    if(sumF2 < SMALL)
                    {
                        continue;
                    }

                    const scalar corr = (1.0 - sumF1)/sumF2;
                    const scalar magSfi = coarseMagSf[i];
                    for (label j=i+1; j<totalNCoarseFaces_; j++)
                    {
                        Fmatrix_()[i][j] *= corr;
                        Fmatrix_()[j][i] = Fmatrix_()[i][j]*magSfi/coarseMagSf[j];
                    }
                }

                {
                    scalar maxErr = -GREAT;
                    scalar minSumF = GREAT;
                    scalar maxSumF = -GREAT;
                    for (label i=0; i<totalNCoarseFaces_; i++)
                    {
                        scalar sumF = 0.0;
                        for (label j=0; j<totalNCoarseFaces_; j++)
                        {
                            sumF += Fmatrix_()[i][j];
                        }
                        maxErr = max(maxErr, mag(sumF - 1.0));
                        minSumF = min(minSumF, sumF);
                        maxSumF = max(maxSumF, sumF);
                    }

                    if(debug)
                    {
                        Info << "max err = " << maxErr
                            << " min sumF = " << minSumF
                            << " max sumF = " << maxSumF << endl;
                    }

                    if(maxErr < 1e-6 || iter == maxIter)
                    {
                        break;
                    }
                }

                iter++;
                for (label i=totalNCoarseFaces_-1; i>=0; i--)
                {
                    scalar sumF1 = 0.0;
                    for (label j=0; j<i; j++)
                    {
                        sumF1 += Fmatrix_()[i][j];
                    }

                    if(sumF1 < SMALL)
                    {
                        continue;
                    }

                    scalar sumF2 = 0.0;
                    for (label j=i+1; j<totalNCoarseFaces_; j++)
                    {
                        sumF2 += Fmatrix_()[i][j];
                    }

                    if(sumF2 > 1.0)
                    {
                        continue;
                    }

                    const scalar corr = (1.0 - sumF2)/sumF1;
                    const scalar magSfi = coarseMagSf[i];
                    for (label j=0; j<i; j++)
                    {
                        Fmatrix_()[i][j] *= corr;
                        Fmatrix_()[j][i] = Fmatrix_()[i][j]*magSfi/coarseMagSf[j];
                    }
                }

                {
                    scalar maxErr = -GREAT;
                    scalar minSumF = GREAT;
                    scalar maxSumF = -GREAT;
                    for (label i=0; i<totalNCoarseFaces_; i++)
                    {
                        scalar sumF = 0.0;
                        for (label j=0; j<totalNCoarseFaces_; j++)
                        {
                            sumF += Fmatrix_()[i][j];
                        }
                        maxErr = max(maxErr, mag(sumF - 1.0));
                        minSumF = min(minSumF, sumF);
                        maxSumF = max(maxSumF, sumF);
                    }

                    if(debug)
                    {
                        Info << "max err = " << maxErr
                            << " min sumF = " << minSumF
                            << " max sumF = " << maxSumF << endl;
                    }

                    if(maxErr < 1e-6 || iter == maxIter)
                    {
                        break;
                    }
                }
            }

            if(iter != maxIter)
            {
                Info << "converged in " << iter << " iterations" << endl;
            }
            else
            {
                Info << "not converged in " << iter << " iterations" << endl;
            }
        }

        constEmissivity_ = readBool(coeffs_.lookup("constantEmissivity"));
        if (constEmissivity_)
        {
            CLU_.reset
            (
                new scalarSquareMatrix
                (
                    totalNCoarseFaces_,
                    0.0
                )
            );

            pivotIndices_.setSize(CLU_().n());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::viewFactor::viewFactor(const volScalarField& T)
:
    radiationModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


Foam::radiation::viewFactor::viewFactor
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    map_(),
    coarseMesh_
    (
        IOobject
        (
            mesh_.name(),
            mesh_.polyMesh::instance(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        finalAgglom_
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fmatrix_(),
    CLU_(),
    selectedPatches_(mesh_.boundary().size(), -1),
    totalNCoarseFaces_(0),
    nLocalCoarseFaces_(0),
    constEmissivity_(false),
    iterCounter_(0),
    pivotIndices_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::viewFactor::~viewFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::viewFactor::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::viewFactor::insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label procI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& Fmatrix
)
{
    forAll (viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(procI, faceI);
        forAll (globalFaces, i)
        {
            Fmatrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}


void Foam::radiation::viewFactor::calculate()
{
    // Store previous iteration
    Qr_.storePrevIter();

    scalarField compactCoarseT(map_->constructSize(), 0.0);
    scalarField compactCoarseE(map_->constructSize(), 0.0);
    scalarField compactCoarseHo(map_->constructSize(), 0.0);

    globalIndex globalNumbering(nLocalCoarseFaces_);

    // Fill local averaged(T), emissivity(E) and external heatFlux(Ho)
    DynamicList<scalar> localCoarseTave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseEave(nLocalCoarseFaces_);
    DynamicList<scalar> localCoarseHoave(nLocalCoarseFaces_);

    forAll (selectedPatches_, i)
    {
        label patchID = selectedPatches_[i];

        const scalarField& Tp = T_.boundaryField()[patchID];
        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        fvPatchScalarField& QrPatch = Qr_.boundaryField()[patchID];

        greyDiffusiveViewFactorFixedValueFvPatchScalarField& Qrp =
            refCast
            <
                greyDiffusiveViewFactorFixedValueFvPatchScalarField
            >(QrPatch);

        const scalarList eb = Qrp.emissivity();

        const scalarList& Hoi = Qrp.Qro();

        const polyPatch& pp = coarseMesh_.boundaryMesh()[patchID];
        const labelList& coarsePatchFace = coarseMesh_.patchFaceMap()[patchID];

        scalarList Tave(pp.size(), 0.0);
        scalarList Eave(Tave.size(), 0.0);
        scalarList Hoiave(Tave.size(), 0.0);

        if (pp.size() > 0)
        {
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom) + 1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            forAll (coarseToFine, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                UIndirectList<scalar> fineSf
                (
                    sf,
                    fineFaces
                );
                scalar area = sum(fineSf());
                // Temperature, emissivity and external flux area weighting
                forAll (fineFaces, j)
                {
                    label faceI = fineFaces[j];
                    Tave[coarseI] += (Tp[faceI]*sf[faceI])/area;
                    Eave[coarseI] += (eb[faceI]*sf[faceI])/area;
                    Hoiave[coarseI] += (Hoi[faceI]*sf[faceI])/area;
                }
            }
        }

        localCoarseTave.append(Tave);
        localCoarseEave.append(Eave);
        localCoarseHoave.append(Hoiave);
    }

    // Fill the local values to distribute
    SubList<scalar>(compactCoarseT, nLocalCoarseFaces_).assign(localCoarseTave);
    SubList<scalar>(compactCoarseE, nLocalCoarseFaces_).assign(localCoarseEave);
    SubList<scalar>
        (compactCoarseHo, nLocalCoarseFaces_).assign(localCoarseHoave);

    // Distribute data
    map_->distribute(compactCoarseT);
    map_->distribute(compactCoarseE);
    map_->distribute(compactCoarseHo);

    // Distribute local global ID
    labelList compactGlobalIds(map_->constructSize(), 0.0);

    labelList localGlobalIds(nLocalCoarseFaces_);

    for(label k = 0; k < nLocalCoarseFaces_; k++)
    {
        localGlobalIds[k] = globalNumbering.toGlobal(Pstream::myProcNo(), k);
    }

    SubList<label>
    (
        compactGlobalIds,
        nLocalCoarseFaces_
    ).assign(localGlobalIds);

    map_->distribute(compactGlobalIds);

    // Create global size vectors
    scalarField T(totalNCoarseFaces_, 0.0);
    scalarField E(totalNCoarseFaces_, 0.0);
    scalarField QrExt(totalNCoarseFaces_, 0.0);

    // Fill lists from compact to global indexes.
    forAll (compactCoarseT, i)
    {
        T[compactGlobalIds[i]] = compactCoarseT[i];
        E[compactGlobalIds[i]] = compactCoarseE[i];
        QrExt[compactGlobalIds[i]] = compactCoarseHo[i];
    }

    Pstream::listCombineGather(T, maxEqOp<scalar>());
    Pstream::listCombineGather(E, maxEqOp<scalar>());
    Pstream::listCombineGather(QrExt, maxEqOp<scalar>());

    Pstream::listCombineScatter(T);
    Pstream::listCombineScatter(E);
    Pstream::listCombineScatter(QrExt);

    // Net radiation
    scalarField q(totalNCoarseFaces_, 0.0);

    if (Pstream::master())
    {
        // Variable emissivity
        if (!constEmissivity_)
        {
            scalarSquareMatrix C(totalNCoarseFaces_, 0.0);

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    scalar invEj = 1.0/E[j];
                    scalar sigmaT4 =
                        radiation::sigmaSB.value()*pow(T[j], 4.0);

                    if (i==j)
                    {
                        C[i][j] = invEj - (invEj - 1.0)*Fmatrix_()[i][j];
                        q[i] += (Fmatrix_()[i][j] - 1.0)*sigmaT4 - QrExt[j];
                    }
                    else
                    {
                        C[i][j] = (1.0 - invEj)*Fmatrix_()[i][j];
                        q[i] += Fmatrix_()[i][j]*sigmaT4 - QrExt[j];
                    }

                }
            }

            Info<< "\nSolving view factor equations (" << CLU_().n()
                << " ..." << flush;
            // Negative coming into the fluid
            scalarSquareMatrix::LUsolve(C, q);
            Info<< "Done." << endl;
        }
        else //Constant emissivity
        {
            // Initial iter calculates CLU and chaches it
            if (iterCounter_ == 0)
            {
                for (label i=0; i<totalNCoarseFaces_; i++)
                {
                    for (label j=0; j<totalNCoarseFaces_; j++)
                    {
                        scalar invEj = 1.0/E[j];
                        if (i==j)
                        {
                            CLU_()[i][j] = invEj-(invEj-1.0)*Fmatrix_()[i][j];
                        }
                        else
                        {
                            CLU_()[i][j] = (1.0 - invEj)*Fmatrix_()[i][j];
                        }
                    }
                }
                Info<< "\nDecomposing C matrix of size " << CLU_().n()
                    << " ... " << endl;
                scalarSquareMatrix::LUDecompose(CLU_(), pivotIndices_);
                Info<< "Done." << endl;
            }

            for (label i=0; i<totalNCoarseFaces_; i++)
            {
                for (label j=0; j<totalNCoarseFaces_; j++)
                {
                    scalar sigmaT4 =
                        radiation::sigmaSB.value()
                       *pow(T[j], 4.0);

                    if (i==j)
                    {
                        q[i] += (Fmatrix_()[i][j] - 1.0)*sigmaT4 - QrExt[j];
                    }
                    else
                    {
                        q[i] += Fmatrix_()[i][j]*sigmaT4 - QrExt[j];
                    }
                }
            }

            Info<< "\nLU Back substitute C matrix... " << endl;
            scalarSquareMatrix::LUBacksubstitute(CLU_(), pivotIndices_, q);
            Info<< "Done." << endl;
            iterCounter_ ++;
        }
    }

    // Scatter q and fill Qr
    Pstream::listCombineScatter(q);
    Pstream::listCombineGather(q, maxEqOp<scalar>());


    label globCoarseId = 0;
    forAll (selectedPatches_, i)
    {
        const label patchID = selectedPatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];
        if (pp.size() > 0)
        {
            scalarField& Qrp = Qr_.boundaryField()[patchID];
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            const labelList& agglom = finalAgglom_[patchID];
            label nAgglom = max(agglom)+1;

            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));

            const labelList& coarsePatchFace =
                coarseMesh_.patchFaceMap()[patchID];

            scalar heatFlux = 0.0;
            forAll (coarseToFine, coarseI)
            {
                label globalCoarse =
                    globalNumbering.toGlobal(Pstream::myProcNo(), globCoarseId);
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceID];
                forAll (fineFaces, k)
                {
                    label faceI = fineFaces[k];

                    Qrp[faceI] = q[globalCoarse];
                    heatFlux += Qrp[faceI]*sf[faceI];
                }
                globCoarseId ++;
            }
        }
    }

    if (debug)
    {
        forAll (Qr_.boundaryField(), patchID)
        {
            const scalarField& Qrp = Qr_.boundaryField()[patchID];
            const scalarField& magSf = mesh_.magSf().boundaryField()[patchID];
            scalar heatFlux = gSum(Qrp*magSf);
            Info<< "Total heat flux at patch: "
                << patchID << " "
                << heatFlux << " [W]" << endl;
        }
    }

    // Relax Qr if necessary
    Qr_.relax();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::viewFactor::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::viewFactor::Ru() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );
}

// ************************************************************************* //
