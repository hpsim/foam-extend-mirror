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

#include "hhuMixtureThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hhuMixtureThermo<MixtureType>::hhuMixtureThermo
(
    const fvMesh& mesh,
    const objectRegistry& obj
)
:
    hhuCombustionThermo(mesh, obj),
    MixtureType(*this, mesh, obj)
{
    scalarField& hCells = h_.internalField();
    scalarField& huCells = hu_.internalField();
    const scalarField& TCells = T_.internalField();
    const scalarField& TuCells = Tu_.internalField();

    forAll(hCells, celli)
    {
        hCells[celli] = this->cellMixture(celli).H(TCells[celli]);
        huCells[celli] = this->cellReactants(celli).H(TuCells[celli]);
    }

    forAll(h_.boundaryField(), patchi)
    {
        h_.boundaryField()[patchi] == h(T_.boundaryField()[patchi], patchi);

        fvPatchScalarField& phu = hu_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];

        forAll(phu, facei)
        {
            phu[facei] = this->patchFaceReactants(patchi, facei).H(pTu[facei]);
        }
    }

    hBoundaryCorrection(h_);
    huBoundaryCorrection(hu_);

    calculate();
    psi_.oldTime();   // Switch on saving old time
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hhuMixtureThermo<MixtureType>::~hhuMixtureThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::hhuMixtureThermo<MixtureType>::calculate()
{
    const scalarField& hCells = h_.internalField();
    const scalarField& huCells = hu_.internalField();
    const scalarField& pCells = p_.internalField();

    scalarField& TCells = T_.internalField();
    scalarField& TuCells = Tu_.internalField();
    scalarField& psiCells = psi_.internalField();
    scalarField& muCells = mu_.internalField();
    scalarField& alphaCells = alpha_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.TH(hCells[celli], TCells[celli]);
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(TCells[celli]);
        alphaCells[celli] = mixture_.alpha(TCells[celli]);

        TuCells[celli] =
            this->cellReactants(celli).TH(huCells[celli], TuCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = p_.boundaryField()[patchi];
        fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = psi_.boundaryField()[patchi];

        fvPatchScalarField& ph = h_.boundaryField()[patchi];
        fvPatchScalarField& phu = hu_.boundaryField()[patchi];

        fvPatchScalarField& pmu_ = mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha_ = alpha_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = mixture_.H(pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture_.mu(pT[facei]);
                palpha_[facei] = mixture_.alpha(pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.TH(ph[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu_[facei] = mixture_.mu(pT[facei]);
                palpha_[facei] = mixture_.alpha(pT[facei]);

                pTu[facei] =
                    this->patchFaceReactants(patchi, facei)
                   .TH(phu[facei], pTu[facei]);
            }
        }
    }
}


template<class MixtureType>
void Foam::hhuMixtureThermo<MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering hhuMixtureThermo<MixtureType>::correct()" << endl;
    }

    // force the saving of the old-time values
    psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "exiting hhuMixtureThermo<MixtureType>::correct()" << endl;
    }
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::hc() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            h_.dimensions()
        )
    );

    volScalarField& hcf = thc();
    scalarField& hcCells = hcf.internalField();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    forAll(hcf.boundaryField(), patchi)
    {
        scalarField& hcp = hcf.boundaryField()[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hhuMixtureThermo<MixtureType>::h
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, celli)
    {
        h[celli] = this->cellMixture(cells[celli]).H(T[celli]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hhuMixtureThermo<MixtureType>::h
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> th(new scalarField(T.size()));
    scalarField& h = th();

    forAll(T, facei)
    {
        h[facei] = this->patchFaceMixture(patchi, facei).H(T[facei]);
    }

    return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hhuMixtureThermo<MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));

    scalarField& cp = tCp();

    forAll(T, facei)
    {
        cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hhuMixtureThermo<MixtureType>::Cp
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp();

    forAll(T, celli)
    {
        cp[celli] = this->cellMixture(cells[celli]).Cp(T[celli]);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::Cp() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cp = tCp();
    scalarField& cpCells = cp.internalField();
    const scalarField& TCells = T_.internalField();

    forAll(TCells, celli)
    {
        cpCells[celli] = this->cellMixture(celli).Cp(TCells[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
    }

    return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hhuMixtureThermo<MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));

    scalarField& cv = tCv();

    forAll(T, facei)
    {
        cv[facei] = this->patchFaceMixture(patchi, facei).Cv(T[facei]);
    }

    return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hhuMixtureThermo<MixtureType>::Cv() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -2, -1, 0)
        )
    );

    volScalarField& cv = tCv();

    forAll(T_, celli)
    {
        cv[celli] = this->cellMixture(celli).Cv(T_[celli]);
    }

    forAll(T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] = Cv(T_.boundaryField()[patchi], patchi);
    }

    return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hhuMixtureThermo<MixtureType>::hu
(
    const scalarField& Tu,
    const labelList& cells
) const
{
    tmp<scalarField> thu(new scalarField(Tu.size()));
    scalarField& hu = thu();

    forAll(Tu, celli)
    {
        hu[celli] = this->cellReactants(cells[celli]).H(Tu[celli]);
    }

    return thu;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::hhuMixtureThermo<MixtureType>::hu
(
    const scalarField& Tu,
    const label patchi
) const
{
    tmp<scalarField> thu(new scalarField(Tu.size()));
    scalarField& hu = thu();

    forAll(Tu, facei)
    {
        hu[facei] = this->patchFaceReactants(patchi, facei).H(Tu[facei]);
    }

    return thu;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::Tb() const
{
    tmp<volScalarField> tTb
    (
        new volScalarField
        (
            IOobject
            (
                "Tb",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_
        )
    );

    volScalarField& Tb_ = tTb();
    scalarField& TbCells = Tb_.internalField();
    const scalarField& TCells = T_.internalField();
    const scalarField& hCells = h_.internalField();

    forAll(TbCells, celli)
    {
        TbCells[celli] =
            this->cellProducts(celli).TH(hCells[celli], TCells[celli]);
    }

    forAll(Tb_.boundaryField(), patchi)
    {
        fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        const fvPatchScalarField& ph = h_.boundaryField()[patchi];
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];

        forAll(pTb, facei)
        {
            pTb[facei] =
                this->patchFaceProducts(patchi, facei)
               .TH(ph[facei], pT[facei]);
        }
    }

    return tTb;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::psiu() const
{
    tmp<volScalarField> tpsiu
    (
        new volScalarField
        (
            IOobject
            (
                "psiu",
                psi_.time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            psi_.dimensions()
        )
    );

    volScalarField& psiu = tpsiu();
    scalarField& psiuCells = psiu.internalField();
    const scalarField& TuCells = Tu_.internalField();
    const scalarField& pCells = p_.internalField();

    forAll(psiuCells, celli)
    {
        psiuCells[celli] =
            this->cellReactants(celli).psi(pCells[celli], TuCells[celli]);
    }

    forAll(psiu.boundaryField(), patchi)
    {
        fvPatchScalarField& ppsiu = psiu.boundaryField()[patchi];

        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];

        forAll(ppsiu, facei)
        {
            ppsiu[facei] =
                this->
                patchFaceReactants(patchi, facei).psi(pp[facei], pTu[facei]);
        }
    }

    return tpsiu;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::psib() const
{
    tmp<volScalarField> tpsib
    (
        new volScalarField
        (
            IOobject
            (
                "psib",
                psi_.time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            psi_.dimensions()
        )
    );

    volScalarField& psib = tpsib();
    scalarField& psibCells = psib.internalField();
    volScalarField Tb_ = Tb();
    const scalarField& TbCells = Tb_.internalField();
    const scalarField& pCells = p_.internalField();

    forAll(psibCells, celli)
    {
        psibCells[celli] =
            this->cellReactants(celli).psi(pCells[celli], TbCells[celli]);
    }

    forAll(psib.boundaryField(), patchi)
    {
        fvPatchScalarField& ppsib = psib.boundaryField()[patchi];

        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        forAll(ppsib, facei)
        {
            ppsib[facei] =
                this->patchFaceReactants
                (patchi, facei).psi(pp[facei], pTb[facei]);
        }
    }

    return tpsib;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::muu() const
{
    tmp<volScalarField> tmuu
    (
        new volScalarField
        (
            IOobject
            (
                "muu",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& muu_ = tmuu();
    scalarField& muuCells = muu_.internalField();
    const scalarField& TuCells = Tu_.internalField();

    forAll(muuCells, celli)
    {
        muuCells[celli] = this->cellReactants(celli).mu(TuCells[celli]);
    }

    forAll(muu_.boundaryField(), patchi)
    {
        fvPatchScalarField& pMuu = muu_.boundaryField()[patchi];
        const fvPatchScalarField& pTu = Tu_.boundaryField()[patchi];

        forAll(pMuu, facei)
        {
            pMuu[facei] =
                this->patchFaceReactants(patchi, facei).mu(pTu[facei]);
        }
    }

    return tmuu;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::hhuMixtureThermo<MixtureType>::mub() const
{
    tmp<volScalarField> tmub
    (
        new volScalarField
        (
            IOobject
            (
                "mub",
                T_.time().timeName(),
                T_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_.mesh(),
            dimensionSet(1, -1, -1, 0, 0)
        )
    );

    volScalarField& mub_ = tmub();
    scalarField& mubCells = mub_.internalField();
    volScalarField Tb_ = Tb();
    const scalarField& TbCells = Tb_.internalField();

    forAll(mubCells, celli)
    {
        mubCells[celli] = this->cellProducts(celli).mu(TbCells[celli]);
    }

    forAll(mub_.boundaryField(), patchi)
    {
        fvPatchScalarField& pMub = mub_.boundaryField()[patchi];
        const fvPatchScalarField& pTb = Tb_.boundaryField()[patchi];

        forAll(pMub, facei)
        {
            pMub[facei] =
                this->patchFaceProducts(patchi, facei).mu(pTb[facei]);
        }
    }

    return tmub;
}


template<class MixtureType>
bool Foam::hhuMixtureThermo<MixtureType>::read()
{
    if (hhuCombustionThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
