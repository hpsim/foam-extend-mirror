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

InClass
    Foam::basicMixture

Description

\*---------------------------------------------------------------------------*/

#ifndef makeBasicMixture_H
#define makeBasicMixture_H

#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBasicMixture(Mixture,Transport,Thermo,EqnOfState)                 \
                                                                              \
typedef Mixture<Transport<specieThermo<Thermo<EqnOfState> > > >               \
    Mixture##Transport##Thermo##EqnOfState;                                   \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
    (Mixture##Transport##Thermo##EqnOfState,                                  \
    #Mixture"<"#Transport"<specieThermo<"#Thermo"<"#EqnOfState">>>>", 0);


#define makeBasicMixturePhys(Mixture,ThermoPhys)                              \
                                                                              \
typedef Mixture<ThermoPhys>   Mixture##ThermoPhys;                            \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Mixture##ThermoPhys,                                                      \
    #Mixture"<"#ThermoPhys">",                                                \
    0                                                                         \
);

#define makeBasicRealFluidMixture(Mixture,Transport,SpecieThermo,Thermo,EqnOfState) \
                                                                              \
typedef Mixture<Transport<SpecieThermo<Thermo<EqnOfState> > > >               \
    Mixture##Transport##SpecieThermo##Thermo##EqnOfState;                     \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
    (Mixture##Transport##SpecieThermo##Thermo##EqnOfState,                    \
    #Mixture"<"#Transport"<"#SpecieThermo"<"#Thermo"<"#EqnOfState">>>>", 0)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBasicRealFluidMixture(Mixture,Transport,SpecieThermo,Thermo,EqnOfState)                 \
                                                                              \
typedef Mixture<Transport<SpecieThermo<Thermo<EqnOfState> > > >               \
    Mixture##Transport##SpecieThermo##Thermo##EqnOfState;                     \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
    (Mixture##Transport##SpecieThermo##Thermo##EqnOfState,                    \
    #Mixture"<"#Transport"<"#SpecieThermo"<"#Thermo"<"#EqnOfState">>>>", 0)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
