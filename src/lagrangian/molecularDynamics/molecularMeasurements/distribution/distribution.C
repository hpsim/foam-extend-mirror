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

\*----------------------------------------------------------------------------*/

#include "distribution.H"
#include "labelList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distribution::distribution()
:
    Map<label>(),
    binWidth_(1)
{}


distribution::distribution(const scalar binWidth)
:
    Map<label>(),
    binWidth_(binWidth)
{}


distribution::distribution(const distribution& d)
:
    Map<label>(static_cast< Map<label> >(d)),
    binWidth_(d.binWidth())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

distribution::~distribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label distribution::totalEntries() const
{
    label sumOfEntries = 0;

    forAllConstIter(Map<label>, *this, iter)
    {
        sumOfEntries += iter();

        if (sumOfEntries < 0)
        {
            WarningIn("label distribution::totalEntries()")
                << "Accumulated distribution values total has become negative: "
                << "sumOfEntries = " << sumOfEntries
                << ". This is most likely to be because too many samples "
                << "have been added to the bins and the label has 'rolled "
                << "round'. Try distribution::approxTotalEntries which "
                << "returns a scalar." << endl;

            sumOfEntries = -1;

            break;
        }
    }

    return sumOfEntries;
}


scalar distribution::approxTotalEntries() const
{
    scalar sumOfEntries = 0;

    forAllConstIter(Map<label>, *this, iter)
    {
        sumOfEntries += scalar(iter());
    }

    return sumOfEntries;
}


scalar distribution::mean() const
{
    scalar runningSum = 0;

    scalar totEnt = approxTotalEntries();

    labelList keys = toc();

    forAll(keys,k)
    {
        label key = keys[k];

        runningSum +=
            (0.5 + scalar(key))
           *binWidth_
           *scalar((*this)[key])
           /totEnt;
    }

    return runningSum;
}


scalar distribution::median()
{
    // From:
    // http://mathworld.wolfram.com/StatisticalMedian.html
    // The statistical median is the value of the distribution variable
    // where the cumulative distribution = 0.5.

    scalar median = 0.0;

    scalar runningSum = 0.0;

    List<Pair<scalar> > normDist(normalised());

    if (normDist.size())
    {
        if (normDist.size() == 1)
        {
            median = normDist[0].first();
        }
        else if
        (
            normDist.size() > 1
         && normDist[0].second()*binWidth_ > 0.5
        )
        {
            scalar xk = normDist[1].first();
            scalar xkm1 = normDist[0].first();
            scalar Sk =
                (normDist[0].second() + normDist[1].second())*binWidth_;
            scalar Skm1 = normDist[0].second()*binWidth_;

            median = (0.5 - Skm1)*(xk - xkm1)/(Sk - Skm1) + xkm1;
        }
        else
        {
            label lastNonZeroIndex = 0;

            forAll(normDist,nD)
            {
                if (runningSum + (normDist[nD].second()*binWidth_) > 0.5)
                {
                    scalar xk = normDist[nD].first();
                    scalar xkm1 = normDist[lastNonZeroIndex].first();
                    scalar Sk = runningSum + (normDist[nD].second()*binWidth_);
                    scalar Skm1 = runningSum;

                    median = (0.5 - Skm1)*(xk - xkm1)/(Sk - Skm1) + xkm1;

                    break;
                }
                else if (normDist[nD].second() > 0.0)
                {
                    runningSum += normDist[nD].second()*binWidth_;

                    lastNonZeroIndex = nD;
                }
            }
        }
    }

    return median;
}


void distribution::add(const scalar valueToAdd)
{
    iterator iter(this->begin());

    label n = label(valueToAdd/binWidth_) - label(neg(valueToAdd/binWidth_));

    iter = find(n);

    if (iter == this->end())
    {
        this->insert(n,1);
    }
    else
    {
        (*this)[n]++;
    }

    if ((*this)[n] < 0)
    {
        FatalErrorIn("distribution::add(const scalar valueToAdd)")
            << "Accumulated distribution value has become negative: "
            << "bin = " << (0.5 + scalar(n)) * binWidth_
            << ", value = " << (*this)[n]
            << ". This is most likely to be because too many samples "
            << "have been added to a bin and the label has 'rolled round'"
            << abort(FatalError);
    }
}


void distribution::add(const label valueToAdd)
{
    add(scalar(valueToAdd));
}


void distribution::insertMissingKeys()
{
    iterator iter(this->begin());

    labelList keys = toc();

    sort(keys);

    if (keys.size())
    {
        for (label k = keys[0]; k < keys[keys.size()-1]; k++)
        {
            iter = find(k);

            if (iter == this->end())
            {
                this->insert(k,0);
            }
        }
    }
}


List< Pair<scalar> > distribution::normalised()
{
    scalar totEnt = approxTotalEntries();

    insertMissingKeys();

    labelList keys = toc();

    sort(keys);

    List<Pair<scalar> > normDist(size());

    forAll(keys,k)
    {
        label key = keys[k];

        normDist[k].first() = (0.5 + scalar(key))*binWidth_;

        normDist[k].second() = scalar((*this)[key])/totEnt/binWidth_;
    }

    return normDist;
}


List< Pair<scalar> > distribution::normalisedMinusMean()
{
    return normalisedShifted(mean());
}


List< Pair<scalar> > distribution::normalisedShifted(const scalar shiftValue)
{
    List<Pair<scalar> > oldDist(normalised());

    List<Pair<scalar> > newDist(oldDist.size());

    forAll(oldDist,u)
    {
        oldDist[u].first() -= shiftValue;
    }

    scalar lowestOldBin = oldDist[0].first()/binWidth_ - 0.5;

    label lowestNewKey = label
    (
        lowestOldBin + 0.5*sign(lowestOldBin)
    );

    scalar interpolationStartDirection =
        sign(scalar(lowestNewKey) - lowestOldBin);

    label newKey = lowestNewKey;

//     Info << shiftValue
//         << nl << lowestOldBin
//         << nl << lowestNewKey
//         << nl << interpolationStartDirection
//         << endl;

//     scalar checkNormalisation = 0;

//     forAll (oldDist, oD)
//     {
//         checkNormalisation += oldDist[oD].second()*binWidth_;
//     }

//     Info << "Initial normalisation = " << checkNormalisation << endl;

    forAll(oldDist,u)
    {
        newDist[u].first() = (0.5 + scalar(newKey)) * binWidth_;

        if (interpolationStartDirection < 0)
        {
            if (u == 0)
            {
                newDist[u].second() =
                    (0.5 + scalar(newKey))*oldDist[u].second()
                  - oldDist[u].second()
                        *(oldDist[u].first() - binWidth_)/ binWidth_;
            }
            else
            {
                newDist[u].second() =
                    (0.5 + scalar(newKey))
                   *(oldDist[u].second() - oldDist[u-1].second())
                  +
                    (
                        oldDist[u-1].second()*oldDist[u].first()
                      - oldDist[u].second()*oldDist[u-1].first()
                    )
                    /binWidth_;
            }
        }
        else
        {
            if (u == oldDist.size() - 1)
            {
                newDist[u].second() =
                    (0.5 + scalar(newKey))*-oldDist[u].second()
                  + oldDist[u].second()*(oldDist[u].first() + binWidth_)
                   /binWidth_;
            }
            else
            {
                newDist[u].second() =
                    (0.5 + scalar(newKey))
                   *(oldDist[u+1].second() - oldDist[u].second())
                  +
                    (
                        oldDist[u].second()*oldDist[u+1].first()
                      - oldDist[u+1].second()*oldDist[u].first()
                    )
                   /binWidth_;
            }
        }

        newKey++;
    }

//     checkNormalisation = 0;

//     forAll (newDist, nD)
//     {
//         checkNormalisation += newDist[nD].second()*binWidth_;
//     }

//     Info << "Shifted normalisation = " << checkNormalisation << endl;

    return newDist;
}


List<Pair<scalar> > distribution::raw()
{
    insertMissingKeys();

    labelList keys = toc();

    sort(keys);

    List<Pair<scalar> > rawDist(size());

    forAll(keys,k)
    {
        label key = keys[k];

        rawDist[k].first() = (0.5 + scalar(key))*binWidth_;

        rawDist[k].second() = scalar((*this)[key]);
    }

    return rawDist;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void distribution::operator=(const distribution& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("distribution::operator=(const distribution&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    Map<label>::operator=(rhs);

    binWidth_ = rhs.binWidth();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const distribution& d)
{
    os  << d.binWidth_
        << static_cast<const Map<label>&>(d);

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const distribution&)"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
