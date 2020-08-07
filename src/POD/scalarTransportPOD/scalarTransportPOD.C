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
    scalarTransportPOD

\*---------------------------------------------------------------------------*/

#include "scalarTransportPOD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scalarTransportPOD, 0);

    addToRunTimeSelectionTable
    (
        PODODE,
        scalarTransportPOD,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::scalarTransportPOD::calcOrthoBase() const
{
    if (orthoBasePtr_)
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
    PtrList<volScalarField> fields(Times.size());

    label nSnapshots = 0;

    forAll (Times, i)
    {
        if (Times[i].equal(0))
        {
            Info << "Skipping time " << Times[i] << endl;

            continue;
        }

        runTime.setTime(Times[i], i);

        // Both field and flux need to be read for the snapshot to be valid
        // Check both

        // Field header
        IOobject fieldHeader
        (
            fieldName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        // Flux header
        IOobject phiHeader
        (
            phiName_,
            runTime.timeName(),
            this->mesh(),
            IOobject::MUST_READ
        );

        if (fieldHeader.headerOk() && phiHeader.headerOk())
        {
            // Record time as valid
            validTimes[i] = true;

            if (firstReadTimeIndex == -1)
            {
                firstReadTimeIndex = i;
            }

            Info<< "    Reading " << fieldName_ << " from time "
                << runTime.timeName() << endl;

            fields.set(nSnapshots, new volScalarField(fieldHeader, mesh()));

            fields[nSnapshots].rename(fieldName_ + name(i));

            if (!reconFieldPtr_)
            {
                Info<< "Reading " << "recon" << fieldName_ << endl;
                reconFieldPtr_ =
                    new volScalarField
                    (
                        "recon" + fieldName_,
                        fields[nSnapshots]
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

    fields.setSize(nSnapshots);

    // Create ortho-normal base for transported variable
    orthoBasePtr_ = new scalarPODOrthoNormalBase(fields, accuracy);

    // Check orthogonality and magnitude of snapshots
    orthoBasePtr_->checkBase();

    Info<< "Write reconstructed snapshots: check" << endl;

    // Reset counter
    nSnapshots = 0;

    forAll (validTimes, i)
    {
        if (validTimes[i])
        {
            runTime.setTime(Times[i], i);

            Info<< "Time = " << runTime.timeName() << endl;

            volScalarField reconField
            (
                IOobject
                (
                    "directRecon" + fieldName_,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "zero",
                    orthoBase().orthoField(0).dimensions(),
                    0
                )
            );

            // Note: use raw pointer access to orthoBase,
            // as it has just been calculated
            for (label i = 0; i < orthoBasePtr_->baseSize(); i++)
            {
                reconField +=
                    orthoBasePtr_->interpolationCoeffs()[nSnapshots][i]*
                    orthoBasePtr_->orthoField(i);
            }

            // Internal field is set.  Correct boundary conditions
            reconField.correctBoundaryConditions();
            reconField.write();

            nSnapshots++;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[origTimeIndex], origTimeIndex);
}


void Foam::scalarTransportPOD::calcDerivativeCoeffs() const
{
    if (derivativePtr_ || lagrangeDerPtr_ || lagrangeSrcPtr_)
    {
        FatalErrorInFunction
            << "Derivative matrix already calculated"
            << abort(FatalError);
    }

    // Calculate coefficients for differential equation
    // Get times list
    Time& runTime = const_cast<Time&>(this->mesh().time());

    // Remember time index to restore it
    label origTimeIndex = runTime.timeIndex();

    // Read diffusivity

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            this->mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    // Read first available flux field.  Note: flux is fixed for
    // scalar transport

    instantList Times = runTime.times();

    // Flux field
    autoPtr<surfaceScalarField> phiPtr;

    forAll (Times, i)
    {
        if (Times[i].equal(0))
        {
            Info << "Skipping time " << Times[i].name() << endl;

            continue;
        }

        runTime.setTime(Times[i], i);

        // Field header
        IOobject fieldHeader
        (
            fieldName_,
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ
        );

        // Flux header
        IOobject phiHeader
        (
            phiName_,
            runTime.timeName(),
            this->mesh(),
            IOobject::MUST_READ
        );

        if (fieldHeader.headerOk() && phiHeader.headerOk())
        {
            Info<< "    Reading " << phiName_ << " from time "
                << runTime.timeName() << endl;

            phiPtr.set(new surfaceScalarField(phiHeader, this->mesh()));
            break;
        }
    }

    // Reset time index to initial state
    runTime.setTime(Times[origTimeIndex], origTimeIndex);

    if (!phiPtr.valid())
    {
        FatalErrorInFunction
            << "Cannot find flux field: " << phiName_
            << abort(FatalError);
    }

    const surfaceScalarField& phi = phiPtr();

    // Create derivative matrix

    const scalarPODOrthoNormalBase& b = orthoBase();

    // Derivative
    derivativePtr_ = new scalarSquareMatrix(b.baseSize(), 0.0);
    scalarSquareMatrix& derivative = *derivativePtr_;

    // Lagrange multiplier derivative
    lagrangeDerPtr_ = new scalarField(b.baseSize(), scalar(0));
    scalarField& lagrangeDer = *lagrangeDerPtr_;

    // Lagrange multiplier source
    lagrangeSrcPtr_ = new scalarField(b.baseSize(), scalar(0));
    scalarField& lagrangeSrc = *lagrangeSrcPtr_;

    const volScalarField& field = *reconFieldPtr_;

    // Possible scaling error: reconsider.  HJ, 5/Aug/2020

    for (label i = 0; i < b.baseSize(); i++)
    {
        const volScalarField& snapI = b.orthoField(i);

        for (label j = 0; j < b.baseSize(); j++)
        {
            const volScalarField& snapJ = b.orthoField(j);

            // Calculate derivative by moving equation terms to rhs
            derivative[i][j] =
                POD::projection
                (
                    fvc::laplacian
                    (
                        DT, snapJ,
                        "laplacian(" + DT.name() + "," + fieldName_ + ")"
                    ),
                    snapI
                )
              - POD::projection
                (
                    fvc::div
                    (
                        phi, snapJ,
                        "div(" + phiName_ + "," + fieldName_ + ")"
                    ),
                    snapI
                );

            // Lagrange multiplier is calculated on boundaries where
            // reconU fixes value
            forAll (field.boundaryField(), patchI)
            {
                if (field.boundaryField()[patchI].fixesValue())
                {
                    // Note pre-multiplication by beta and signs
                    // of derivative and source.  Su-Sp treatment
                    lagrangeDer[i] += -beta_*
                        POD::projection
                        (
                            snapJ.boundaryField()[patchI],
                            snapI.boundaryField()[patchI]
                        );

                    lagrangeSrc[i] += beta_*
                        POD::projection
                        (
                            field.boundaryField()[patchI],
                            snapI.boundaryField()[patchI]
                        );
                }
            }
        }
    }

    Info<< "derivative: " << derivative << nl
        << "lagrangeDer: " << lagrangeDer << nl
        << "lagrangeSrc: " << lagrangeSrc << endl;
}


void Foam::scalarTransportPOD::updateFields() const
{
    if (fieldUpdateTimeIndex_ < mesh().time().timeIndex())
    {
        // Field update required.  Record update time index
        fieldUpdateTimeIndex_ = mesh().time().timeIndex();

        if (!reconFieldPtr_)
        {
            FatalErrorInFunction
                << "Reconstructed field not allocated"
                << abort(FatalError);
        }

        volScalarField& field = *reconFieldPtr_;

        const scalarPODOrthoNormalBase& b = orthoBase();

        field = dimensionedScalar("zero", b.orthoField(0).dimensions(), 0);
        Info<< "coeffs: " << coeffs_ << endl;
        forAll (coeffs_, i)
        {
            field += coeffs_[i]*b.orthoField(i);
        }

        field.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::scalarTransportPOD::scalarTransportPOD
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PODODE(mesh, dict),
    fieldName_(dict.lookup("field")),
    phiName_(dict.lookup("flux")),
    beta_(readScalar(dict.lookup("beta"))),
    coeffs_(),
    derivativePtr_(nullptr),
    lagrangeDerPtr_(nullptr),
    lagrangeSrcPtr_(nullptr),
    orthoBasePtr_(nullptr),
    reconFieldPtr_(nullptr),
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

Foam::scalarTransportPOD::~scalarTransportPOD()
{
    deleteDemandDrivenData(derivativePtr_);
    deleteDemandDrivenData(lagrangeDerPtr_);
    deleteDemandDrivenData(lagrangeSrcPtr_);

    deleteDemandDrivenData(orthoBasePtr_);

    deleteDemandDrivenData(reconFieldPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::scalarTransportPOD::nEqns() const
{
    return coeffs().size();
}


Foam::scalarField& Foam::scalarTransportPOD::coeffs()
{
    return coeffs_;
}


const Foam::scalarField& Foam::scalarTransportPOD::coeffs() const
{
    return coeffs_;
}


void Foam::scalarTransportPOD::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    if (!derivativePtr_ || !lagrangeDerPtr_ || !lagrangeSrcPtr_)
    {
        calcDerivativeCoeffs();
    }

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
        }
    }
}


void Foam::scalarTransportPOD::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    derivatives(x, y, dfdx);

    // Clear jacobian matrix
    // Note: this needs to be done component-by-components
    forAll (y, i)
    {
        forAll (y, j)
        {
            dfdy[i][j] = 0;
        }
    }
}


const Foam::scalarPODOrthoNormalBase&
Foam::scalarTransportPOD::orthoBase() const
{
    if (!orthoBasePtr_)
    {
        calcOrthoBase();
    }

    return *orthoBasePtr_;
}


const Foam::volScalarField& Foam::scalarTransportPOD::reconField() const
{
    updateFields();

    return *reconFieldPtr_;
}


void Foam::scalarTransportPOD::writeSnapshots() const
{
    const scalarPODOrthoNormalBase& b = orthoBase();

    Info<< "Writing POD base for Time = " << mesh().time().timeName() << endl;

    for (label i = 0; i < b.baseSize(); i++)
    {
        b.orthoField(i).write();
    }
}


void Foam::scalarTransportPOD::write() const
{
    // Recalculate field and force a write
    updateFields();
    reconField().write();
}


// ************************************************************************* //
