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

inline Foam::PtrList<Foam::volScalarField>&
Foam::basicMultiComponentMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basicMultiComponentMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::basicMultiComponentMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::basicMultiComponentMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::basicMultiComponentMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basicMultiComponentMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}


inline bool Foam::basicMultiComponentMixture::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}


inline Foam::scalar Foam::basicMultiComponentMixture::fres
(
    const scalar ft,
    const scalar stoicRatio
) const
{
    return max(ft - (1.0 - ft)/stoicRatio, 0.0);
}


inline Foam::tmp<Foam::volScalarField> Foam::basicMultiComponentMixture::fres
(
    const volScalarField& ft,
    const dimensionedScalar& stoicRatio
) const
{
    return max(ft - (scalar(1) - ft)/stoicRatio.value(), scalar(0));
}


// ************************************************************************* //
