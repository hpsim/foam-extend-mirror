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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class thermo>
inline constTransport<thermo>::constTransport
(
    const thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    thermo(t),
    Mu(mu),
    rPr(1.0/Pr)
{}


// Construct as named copy
template<class thermo>
inline constTransport<thermo>::constTransport
(
    const word& name,
    const constTransport& ct
)
:
    thermo(name, ct),
    Mu(ct.Mu),
    rPr(ct.rPr)
{}


// Construct and return a clone
template<class thermo>
inline autoPtr<constTransport<thermo> > constTransport<thermo>::clone
() const
{
    return autoPtr<constTransport<thermo> >
    (
        new constTransport<thermo>(*this)
    );
}


// Selector from Istream
template<class thermo>
inline autoPtr<constTransport<thermo> > constTransport<thermo>::New
(
    Istream& is
)
{
    return autoPtr<constTransport<thermo> >
    (
        new constTransport<thermo>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Dynamic viscosity [kg/ms]
template<class thermo>
inline scalar constTransport<thermo>::mu(const scalar) const
{
    return Mu;
}


// Thermal conductivity [W/mK]
template<class thermo>
inline scalar constTransport<thermo>::kappa(const scalar T) const
{
    return this->Cp(T)*mu(T)*rPr;
}


// Thermal diffusivity for enthalpy [kg/ms]
template<class thermo>
inline scalar constTransport<thermo>::alpha(const scalar T) const
{
    scalar Cp_ = this->Cp(T);

    scalar deltaT = T - specie::Tstd();
    scalar CpBar =
        (deltaT*(this->H(T) - this->H(specie::Tstd())) + Cp_)/(sqr(deltaT) + 1);

    return Cp_*mu(T)*rPr/CpBar;
}

// CL: for real gas thermo
// Thermal conductivity [W/mK]
template<class thermo>
inline scalar constTransport<thermo>::kappa(const scalar rho, const scalar T) const
{
    return this->Cp(rho, T)*mu(T)*rPr;
}

// CL: for real gas thermo
// Thermal diffusivity for enthalpy [kg/ms]
template<class thermo>
inline scalar constTransport<thermo>::alpha
(
    const scalar rho,
    const scalar T
) const
{
    scalar Cp_ = this->Cp(rho, T);

    scalar deltaT = T - specie::Tstd();
    scalar CpBar =
    (
        deltaT*(this->H(rho, T) - this->H(this->rhostd(), specie::Tstd()))
      + Cp_
    )/(sqr(deltaT) + 1);

    return Cp_*mu(T)*rPr/CpBar;
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::constTransport<thermo>::operator+=
(
    const constTransport<thermo>& st
)
{
    scalar molr1 = this->nMoles();

    thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    Mu = molr1*Mu + molr2*st.Mu;
    rPr = 1.0/(molr1/rPr + molr2/st.rPr);
}


template<class thermo>
inline void Foam::constTransport<thermo>::operator-=
(
    const constTransport<thermo>& st
)
{
    scalar molr1 = this->nMoles();

    thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    Mu = molr1*Mu - molr2*st.Mu;
    rPr = 1.0/(molr1/rPr - molr2/st.rPr);
}


template<class thermo>
inline void Foam::constTransport<thermo>::operator*=
(
    const scalar s
)
{
    thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline constTransport<thermo> operator+
(
    const constTransport<thermo>& ct1,
    const constTransport<thermo>& ct2
)
{
    thermo t
    (
        static_cast<const thermo&>(ct1) + static_cast<const thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constTransport<thermo>
    (
        t,
        molr1*ct1.Mu + molr2*ct2.Mu,
        molr1*ct1.rPr + molr2*ct2.rPr
    );
}


template<class thermo>
inline constTransport<thermo> operator-
(
    const constTransport<thermo>& ct1,
    const constTransport<thermo>& ct2
)
{
    thermo t
    (
        static_cast<const thermo&>(ct1) - static_cast<const thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constTransport<thermo>
    (
        t,
        molr1*ct1.Mu - molr2*ct2.Mu,
        molr1*ct1.rPr - molr2*ct2.rPr
    );
}


template<class thermo>
inline constTransport<thermo> operator*
(
    const scalar s,
    const constTransport<thermo>& ct
)
{
    return constTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.Mu,
        ct.rPr
    );
}


template<class thermo>
inline constTransport<thermo> operator==
(
    const constTransport<thermo>& ct1,
    const constTransport<thermo>& ct2
)
{
    return ct2 - ct1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
