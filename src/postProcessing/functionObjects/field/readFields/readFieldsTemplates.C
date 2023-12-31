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

#include "readFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::readFields::loadField
(
    const word& fieldName,
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& vflds,
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >& sflds
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfType;

    if (obr_.foundObject<vfType>(fieldName))
    {
        if (debug)
        {
            Info<< "readFields : Field " << fieldName << " already in database"
                << endl;
        }
    }
    else if (obr_.foundObject<sfType>(fieldName))
    {
        if (debug)
        {
            Info<< "readFields : Field " << fieldName << " already in database"
                << endl;
        }
    }
    else
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == vfType::typeName
        )
        {
            // store field locally
            Info<< "    Reading " << fieldName << endl;
            label sz = vflds.size();
            vflds.setSize(sz+1);
            vflds.set(sz, new vfType(fieldHeader, mesh));
        }
        else if
        (
            fieldHeader.headerOk()
         && fieldHeader.headerClassName() == sfType::typeName
        )
        {
            // store field locally
            Info<< "    Reading " << fieldName << endl;
            label sz = sflds.size();
            sflds.setSize(sz+1);
            sflds.set(sz, new sfType(fieldHeader, mesh));
        }
    }
}


// ************************************************************************* //
