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

Class
    blockVectorNSAMGInterfaceFields

Description
    Macros for VectorN types for AMG interface fields with block coeffs

Author
    Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#include "BlockSAMGInterfaceField.H"
#include "ProcessorBlockSAMGInterfaceField.H"
#include "GGIBlockSAMGInterfaceField.H"
#include "VectorNFieldTypes.H"
#include "ExpandTensorNField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                     \
                                                                              \
typedef BlockSAMGInterfaceField<type > block##Type##SAMGInterfaceField;       \
defineNamedTemplateTypeNameAndDebug(block##Type##SAMGInterfaceField, 0);      \
defineTemplateRunTimeSelectionTable(block##Type##SAMGInterfaceField, lduInterface); \
                                                                              \
typedef ProcessorBlockSAMGInterfaceField<type > block##Type##ProcessorSAMGInterfaceField;  \
makeBlockSAMGInterfaceField(block##Type##SAMGInterfaceField, block##Type##ProcessorSAMGInterfaceField); \
                                                                              \
typedef GGIBlockSAMGInterfaceField<type > block##Type##GGISAMGInterfaceField; \
makeBlockSAMGInterfaceField(block##Type##SAMGInterfaceField, block##Type##GGISAMGInterfaceField); \

forAllVectorNTypes(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
