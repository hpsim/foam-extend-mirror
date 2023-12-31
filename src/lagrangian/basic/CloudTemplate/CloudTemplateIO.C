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

#include "CloudTemplate.H"
#include "Particle.H"
#include "foamTime.H"
#include "IOPosition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::word Foam::Cloud<ParticleType>::cloudPropertiesName("cloudProperties");


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::readCloudUniformProperties()
{
    IOobject uniformPropsDictHeader
    (
        cloudPropertiesName,
        time().timeName(),
        "uniform"/cloud::prefix/name(),
        db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (uniformPropsDictHeader.headerOk())
    {
        const IOdictionary uniformPropsDict(uniformPropsDictHeader);

        word procName("processor" + Foam::name(Pstream::myProcNo()));
        if (uniformPropsDict.found(procName))
        {
            uniformPropsDict.subDict(procName).lookup("particleCount")
                >> particleCount_;
        }
    }
    else
    {
        particleCount_ = 0;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeCloudUniformProperties() const
{
    IOdictionary uniformPropsDict
    (
        IOobject
        (
            cloudPropertiesName,
            time().timeName(),
            "uniform"/cloud::prefix/name(),
            db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelList np(Pstream::nProcs(), 0);
    np[Pstream::myProcNo()] = particleCount_;

    Pstream::listCombineGather(np, maxEqOp<label>());
    Pstream::listCombineScatter(np);

    forAll(np, i)
    {
        word procName("processor" + Foam::name(i));
        uniformPropsDict.add(procName, dictionary());
        uniformPropsDict.subDict(procName).add("particleCount", np[i]);
    }

    uniformPropsDict.regIOobject::write();
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::initCloud(const bool checkClass)
{
    readCloudUniformProperties();

    IOPosition<ParticleType> ioP(*this);

    if (ioP.headerOk())
    {
        ioP.readData(*this, checkClass);
        ioP.close();

        if (this->size())
        {
            readFields();
        }
    }
    else
    {
        WarningIn("Cloud<ParticleType>::initCloud(const bool checkClass)")
            << "Cannot read particle positions file " << nl
            << "    " << ioP.path() << nl
            << "    assuming the initial cloud contains 0 particles." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const bool checkClass
)
:
    cloud(pMesh),
    polyMesh_(pMesh),
    particleCount_(0)
{
    initCloud(checkClass);
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const bool checkClass
)
:
    cloud(pMesh, cloudName),
    polyMesh_(pMesh),
    particleCount_(0)
{
    initCloud(checkClass);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
Foam::IOobject Foam::Cloud<ParticleType>::fieldIOobject
(
    const word& fieldName,
    const IOobject::readOption r
) const
{
    return IOobject
    (
        fieldName,
        time().timeName(),
        *this,
        r,
        IOobject::NO_WRITE,
        false
    );
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldIOobject
(
    const Cloud<ParticleType>& c,
    const IOField<DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorIn
        (
            "void Cloud<ParticleType>::checkFieldIOobject"
            "(const Cloud<ParticleType>&, const IOField<DataType>&) const"
        )   << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::readFields()
{}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeFields() const
{
    if (this->size())
    {
        const ParticleType& p = *this->first();
        ParticleType::writeFields(p.cloud());
    }
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    writeCloudUniformProperties();

    if (this->size())
    {
        writeFields();
        return cloud::writeObject(fmt, ver, cmp);
    }
    else
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Cloud<ParticleType>& pc)
{
    pc.writeData(os);

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const Cloud<ParticleType>&)");

    return os;
}


// ************************************************************************* //
