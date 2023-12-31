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

#include "IOPosition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::IOPosition<ParticleType>::IOPosition
(
    const Cloud<ParticleType>& c
)
:
    regIOobject
    (
        IOobject
        (
            "positions",
            c.time().timeName(),
            c,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    cloud_(c)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
bool Foam::IOPosition<ParticleType>::write() const
{
    if (cloud_.size())
    {
        return regIOobject::write();
    }
    else
    {
        return true;
    }
}


template<class ParticleType>
bool Foam::IOPosition<ParticleType>::writeData(Ostream& os) const
{
    os<< cloud_.size() << nl << token::BEGIN_LIST << nl;

    forAllConstIter(typename Cloud<ParticleType>, cloud_, iter)
    {
        // Prevent writing additional fields
        static_cast<const Particle<ParticleType>&>(iter()).write
        (
            os,
            false
        );
        os  << nl;
    }

    os<< token::END_LIST << endl;

    return os.good();
}


template<class ParticleType>
void Foam::IOPosition<ParticleType>::readData
(
    Cloud<ParticleType>& c,
    bool checkClass
)
{
    Istream& is = readStream(checkClass ? typeName : "");

    token firstToken(is);

    if (firstToken.isLabel())
    {
        label s = firstToken.labelToken();

        // Read beginning of contents
        is.readBeginList("Cloud<ParticleType>");

        for (label i=0; i<s; i++)
        {
            // Do not read any fields, position only
            c.append(new ParticleType(c, is, false));
        }

        // Read end of contents
        is.readEndList("Cloud<ParticleType>");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "void IOPosition<ParticleType>::readData"
                "(Cloud<ParticleType>&, bool)",
                is
            )   << "incorrect first token, '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);
            // Do not read any fields, position only
            c.append(new ParticleType(c, is, false));
            is >> lastToken;
        }
    }
    else
    {
        FatalIOErrorIn
        (
            "void IOPosition<ParticleType>::readData"
            "(Cloud<ParticleType>&, bool)",
            is
        )   << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    // Check state of IOstream
    is.check
    (
        "void IOPosition<ParticleType>::readData(Cloud<ParticleType>&, bool)"
    );
}


// ************************************************************************* //
