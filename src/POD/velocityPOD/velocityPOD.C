/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Class
    velocityPOD

\*---------------------------------------------------------------------------*/

#include "velocityPOD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityPOD, 0);

    addToRunTimeSelectionTable
    (
        PODODE,
        velocityPOD,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::velocityPOD::calcOrthoBase() const
{
    if
    (
        orthoBasePtr_
     || pBasePtr_
     || phiBasePtr_
     || reconUPtr_
     || reconPPtr_
    )
    {
        FatalErrorInFunction
            << "Orthogonal base already calculated"
            << abort(FatalError);
    }

    // Create ortho-normal base
    scalar accuracy = readScalar(dict().lookup("accuracy"));

    // Get times list
    Time& runTime = const_cast<Time&>(mesh().time());

    // Remember time index to restore it after the scan
    label origTimeIndex = runTime.timeIndex();
    label firstReadTimeIndex = -1;

    instantList Times = runTime.times();

    // Assume no times are valid
    boolList validTimes(Times.size(), false);

    // Create a list of snapshots
    PtrList<volVectorField> UFields(Times.size());

    PtrList<volScalarField> pFields(Times.size());

    PtrList<surfaceScalarField> phiFields(Times.size());

    label nSnapshots = 0;

    forAll (Times, i)
    {
        if (Times[i].equal(0))
        {
            Info << "Skipping time " << Times[i].name() << endl;

            continue;
        }

        runTime.setTime(Times[i], i);

        IOobject UHeader
        (
            UName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        IOobject pHeader
        (
            pName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        IOobject phiHeader
        (
            phiName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        if (UHeader.headerOk() && pHeader.headerOk() && phiHeader.headerOk())
        {
            // Record time as valid
            validTimes[i] = true;

            if (firstReadTimeIndex == -1)
            {
                firstReadTimeIndex = i;
            }

            Info<< "Reading snapshots from time = "
                << runTime.timeName() << endl;

            // Get velocity snapshot
            UFields.set
            (
                nSnapshots,
                new volVectorField(UHeader, this->mesh())
            );

            UFields[nSnapshots].rename(UName_ + name(i));

            // Get pressure snapshot
            pFields.set
            (
                nSnapshots,
                new volScalarField(pHeader, this->mesh())
            );

            pFields[nSnapshots].rename(pName_ + name(i));

            // Get flux snapshot
            phiFields.set
            (
                nSnapshots,
                new surfaceScalarField(phiHeader, this->mesh())
            );

            phiFields[nSnapshots].rename(phiName_ + name(i));

            // Read the first time into reconU and reconP for
            // boundary conditions and general reconstruction settings
            if (!reconUPtr_)
            {
                Info<< "Reading " << "recon" << UName_ << endl;
                reconUPtr_ =
                    new volVectorField
                    (
                        "recon" + UName_,
                        UFields[nSnapshots]
                    );
            }

            if (!reconPPtr_)
            {
                Info<< "Reading " << "recon" << pName_ << endl;
                reconPPtr_ =
                    new volScalarField
                    (
                        "recon" + pName_,
                        pFields[nSnapshots]
                    );
            }

            nSnapshots++;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[firstReadTimeIndex], firstReadTimeIndex);

    // Resize snapshots
    if (nSnapshots < 2)
    {
        FatalErrorInFunction
            << "Insufficient number of snapshots: " << nSnapshots
            << abort(FatalError);
    }

    Info << "Number of snapshots: " << nSnapshots << endl;

    UFields.setSize(nSnapshots);
    pFields.setSize(nSnapshots);
    phiFields.setSize(nSnapshots);

    // Create ortho-normal base for velocity
    orthoBasePtr_ = new vectorPODOrthoNormalBase(UFields, accuracy);

    // Check orthogonality and magnitude of snapshots
    orthoBasePtr_->checkBase();

    // Collect pressure field base
    pBasePtr_ = new PtrList<volScalarField>(orthoBasePtr_->baseSize());
    orthoBasePtr_->calcOrthoBase(pFields, *pBasePtr_);

    // Collect flux field base
    phiBasePtr_ = new PtrList<surfaceScalarField>(orthoBasePtr_->baseSize());
    orthoBasePtr_->calcOrthoBase(phiFields, *phiBasePtr_);


    Info<< "Write reconstructed snapshots: check" << endl;

    // Reset counter
    nSnapshots = 0;

    forAll (validTimes, i)
    {
        if (validTimes[i])
        {
            runTime.setTime(Times[i], i);

            Info<< "Time = " << runTime.timeName() << endl;

            const volVectorField& reconU = *reconUPtr_;
            const volScalarField& reconP = *reconPPtr_;

            volVectorField directReconU
            (
                "directReconU",
                reconU
            );

            volScalarField directReconP
            (
                "directReconP",
                reconP
            );

            directReconU =
                dimensionedVector("zero", reconU.dimensions(), vector::zero);

            directReconP =
                dimensionedScalar("zero", reconP.dimensions(), scalar(0));

            vectorField& UIn = directReconU.internalField();
            scalarField& pIn = directReconP.internalField();

            for (label i = 0; i < orthoBasePtr_->baseSize(); i++)
            {
                // Update velocity
                UIn +=
                    orthoBasePtr_->interpolationCoeffs()[nSnapshots][i]*
                    orthoBasePtr_->orthoField(i);

                // Update pressure
                pIn +=
                    orthoBasePtr_->interpolationCoeffs()[nSnapshots][i]*
                    pBasePtr_->operator[](i);
            }

            // Internal field is set.  Correct boundary conditions
            directReconU.correctBoundaryConditions();
            directReconU.write();

            directReconP.correctBoundaryConditions();
            directReconP.write();

            nSnapshots++;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[origTimeIndex], origTimeIndex);
}


void Foam::velocityPOD::calcDerivativeCoeffs() const
{
    if
    (
        convectionDerivativePtr_
     || derivativePtr_
     || lagrangeDerPtr_
     || lagrangeSrcPtr_
    )
    {
        FatalErrorInFunction
            << "Derivative matrix already calculated"
            << abort(FatalError);
    }

    // Calculate coefficients for differential equation

    // Get bases.  Note: velocity base is in ortho base
    const vectorPODOrthoNormalBase& b = orthoBase();

    // Pressure field base
    const PtrList<volScalarField>& pb = pBase();

    // Flux field base
    const PtrList<surfaceScalarField>& phib = phiBase();

    // Create derivative matrices

    // Convection derivative
    convectionDerivativePtr_ = new PtrList<scalarSquareMatrix>(b.baseSize());
    PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    forAll (convectionDerivative, i)
    {
        convectionDerivative.set
        (
            i,
            new scalarSquareMatrix(b.baseSize(), scalar(0))
        );
    }

    // Diffusion and pressure derivative
    derivativePtr_ = new scalarSquareMatrix(b.baseSize(), scalar(0));
    scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative
    lagrangeDerPtr_ = new scalarField(b.baseSize(), scalar(0));
    scalarField& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source
    lagrangeSrcPtr_ = new scalarField(b.baseSize(), scalar(0));
    scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    // Get solution field four boundary conditions.  Force raw access without
    // checking
    const volVectorField& U = *reconUPtr_;

    // i = equation number (normal basis in dot-product)
    // j = summation over all bases
    // k = second summation for convection term

    // Calculate derivative by moving equation terms to rhs

    // Calculate volume for scaling - THIS IS WRONG!!! HJ, HERE!!!
    const scalar V = 0.2*gSum(mesh().V().field());

    // Derivatives assembly loop
    for (label i = 0; i < b.baseSize(); i++)
    {
        const volVectorField& snapUI = b.orthoField(i);

        for (label j = 0; j < b.baseSize(); j++)
        {
            const volVectorField& snapUJ = b.orthoField(j);

            const volScalarField& snapPJ = pb[j];

            derivative[i][j] =
                POD::projection
                (
                    fvc::laplacian
                    (
                        nu_, snapUJ,
                        "laplacian(" + nu_.name() + "," + UName_ + ")"
                    ),
                    snapUI
                )*V
              - POD::projection(fvc::grad(snapPJ, "grad(p)"), snapUI)*V;

            for (label k = 0; k < b.baseSize(); k++)
            {
                const surfaceScalarField& snapPhiK = phib[k];

                convectionDerivative[i][j][k] =
                   -POD::projection
                    (
                        fvc::div
                        (
                            snapPhiK, snapUJ,
                            "div(" + phiName_ + "," + UName_ + ")"
                        ),
                        snapUI
                    )*V;
            }

            // Lagrange multiplier is calculated on boundaries where
            // reconU fixes value
            forAll (U.boundaryField(), patchI)
            {
                if (U.boundaryField()[patchI].fixesValue())
                {
                    // Note pre-multiplication by beta and signs
                    // of derivative and source.  Su-Sp treatment
                    lagrangeDer[i] += -beta_*
                        POD::projection
                        (
                            snapUJ.boundaryField()[patchI],
                            snapUI.boundaryField()[patchI]
                        );

                    lagrangeSrc[i] += beta_*
                        POD::projection
                        (
                            U.boundaryField()[patchI],
                            snapUI.boundaryField()[patchI]
                        );
                }
            }

        }
    }

    Info<< "convectionDerivative: " << convectionDerivative << nl
        << "derivative: " << derivative << nl
        << "lagrangeDer: " << lagrangeDer << nl
        << "lagrangeSrc: " << lagrangeSrc << endl;
}


void Foam::velocityPOD::updateFields() const
{
    if (fieldUpdateTimeIndex_ < mesh().time().timeIndex())
    {
        // Field update required.  Record update time index
        fieldUpdateTimeIndex_ = mesh().time().timeIndex();

        if (!reconUPtr_ || !reconPPtr_)
        {
            FatalErrorInFunction
                << "Reconstructed field not allocated"
                << abort(FatalError);
        }

        volVectorField& reconU = *reconUPtr_;
        volScalarField& reconP = *reconPPtr_;

        reconU = dimensionedVector("zero", reconU.dimensions(), vector::zero);
        reconP = dimensionedScalar("zero", reconP.dimensions(), scalar(0));

        vectorField& reconUIn = reconU.internalField();
        scalarField& reconPIn = reconP.internalField();

        const vectorPODOrthoNormalBase& b = orthoBase();

        const PtrList<volScalarField>& pB = pBase();

        forAll (coeffs_, i)
        {
            // Update velocity
            reconUIn += coeffs_[i]*b.orthoField(i);

            // Update pressure
            reconPIn += coeffs_[i]*pB[i];
        }

        // Internal field is set.  Correct boundary conditions
        reconU.correctBoundaryConditions();
        reconP.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::velocityPOD::velocityPOD
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PODODE(mesh, dict),
    UName_(dict.lookup("U")),
    pName_(dict.lookup("p")),
    phiName_(dict.lookup("phi")),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            this->mesh().time().constant(),
            this->mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nu_(transportProperties_.lookup("nu")),
    beta_(readScalar(dict.lookup("beta"))),
    coeffs_(),
    convectionDerivativePtr_(nullptr),
    derivativePtr_(nullptr),
    lagrangeDerPtr_(nullptr),
    lagrangeSrcPtr_(nullptr),
    orthoBasePtr_(nullptr),
    pBasePtr_(nullptr),
    phiBasePtr_(nullptr),
    reconUPtr_(nullptr),
    reconPPtr_(nullptr),
    fieldUpdateTimeIndex_(-1)
{
    // Check beta
    if (beta_ < 0)
    {
        FatalErrorInFunction
            << "Negative beta: " << beta_
            << abort(FatalError);
    }

    // Grab coefficients from the first snapshot of the ortho-normal base
    coeffs_.setSize(orthoBase().baseSize());

    const scalarRectangularMatrix& orthoBaseCoeffs =
        orthoBase().interpolationCoeffs();

    forAll (coeffs_, i)
    {
        coeffs_[i] = orthoBaseCoeffs[0][i];
    }
    Info<< "Zero coeffs: " << coeffs_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityPOD::~velocityPOD()
{
    deleteDemandDrivenData(convectionDerivativePtr_);
    deleteDemandDrivenData(derivativePtr_);
    deleteDemandDrivenData(lagrangeDerPtr_);
    deleteDemandDrivenData(lagrangeSrcPtr_);

    deleteDemandDrivenData(orthoBasePtr_);
    deleteDemandDrivenData(pBasePtr_);
    deleteDemandDrivenData(phiBasePtr_);

    deleteDemandDrivenData(reconUPtr_);
    deleteDemandDrivenData(reconPPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::velocityPOD::nEqns() const
{
    return coeffs().size();
}


Foam::scalarField& Foam::velocityPOD::coeffs()
{
    return coeffs_;
}


const Foam::scalarField& Foam::velocityPOD::coeffs() const
{
    return coeffs_;
}


void Foam::velocityPOD::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    if
    (
        !convectionDerivativePtr_
     || !derivativePtr_
     || !lagrangeDerPtr_
     || !lagrangeSrcPtr_
    )
    {
        calcDerivativeCoeffs();
    }

    // Convection derivative
    const PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    // Diffusion and pressure derivative
    const scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative, boundary conditions
    const scalarField& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source, boundary conditions
    const scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    forAll (dydx, i)
    {
        dydx[i] = 0;
        dydx[i] = lagrangeSrc[i];

        forAll (y, j)
        {
            dydx[i] += derivative[i][j]*y[j] + lagrangeDer[i]*y[j];

            forAll (y, k)
            {
                dydx[i] += convectionDerivative[i][j][k]*y[j]*y[k];
            }
        }
    }
}


void Foam::velocityPOD::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    if (!convectionDerivativePtr_ || !derivativePtr_)
    {
        calcDerivativeCoeffs();
    }

    // Convection derivative
    const PtrList<scalarSquareMatrix>& convectionDerivative =
        *convectionDerivativePtr_;

    // Diffusion and pressure derivative
    const scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative, boundary conditions
    const scalarField& lagrangeDer = *lagrangeDerPtr_;

    // Must calculate derivatives
    derivatives(x, y, dfdx);

    // Clear jacobian matrix
    forAll (y, i)
    {
        forAll (y, j)
        {
            dfdy[i][j] = 0;
        }
    }

    forAll (y, i)
    {
        forAll (y, j)
        {
            dfdy[i][j] = derivative[i][j] + lagrangeDer[i];

            forAll (y, k)
            {
                dfdy[i][j] += convectionDerivative[i][j][k]*y[k];
            }
        }
    }
}


void Foam::velocityPOD::update(const scalar delta)
{
    // Warning: lots of cost here: updating fields for top-level
    // function objects.  Can be removed from production code
    // HJ, 5/Aug/2020
    updateFields();
}


const Foam::PODOrthoNormalBase<Foam::vector>&
Foam::velocityPOD::orthoBase() const
{
    if (!orthoBasePtr_)
    {
        calcOrthoBase();
    }

    return *orthoBasePtr_;
}


const Foam::PtrList<Foam::volScalarField>&
Foam::velocityPOD::pBase() const
{
    if (!pBasePtr_)
    {
        calcOrthoBase();
    }

    return *pBasePtr_;
}


const Foam::PtrList<Foam::surfaceScalarField>&
Foam::velocityPOD::phiBase() const
{
    if (!phiBasePtr_)
    {
        calcOrthoBase();
    }

    return *phiBasePtr_;
}


const Foam::volVectorField& Foam::velocityPOD::reconU() const
{
    updateFields();

    return *reconUPtr_;
}


const Foam::volScalarField& Foam::velocityPOD::reconP() const
{
    updateFields();

    return *reconPPtr_;
}


void Foam::velocityPOD::writeSnapshots() const
{
    const vectorPODOrthoNormalBase& b = orthoBase();
    const PtrList<volScalarField>& pb = pBase();

    Info<< "Writing POD base for Time = " << mesh().time().timeName() << endl;

    for (label baseI = 0; baseI < b.baseSize(); baseI++)
    {
        // Write snapshots
        volVectorField snapU
        (
            UName_ + "POD" + name(baseI),
            b.orthoField(baseI)
        );
        snapU.write();

        volScalarField snapP
        (
            pName_ + "POD" + name(baseI),
            pb[baseI]
        );
        snapP.write();
    }
}


void Foam::velocityPOD::write() const
{
    // Recalculate field and force a write
    updateFields();
    reconU().write();
    reconP().write();
}


// ************************************************************************* //
