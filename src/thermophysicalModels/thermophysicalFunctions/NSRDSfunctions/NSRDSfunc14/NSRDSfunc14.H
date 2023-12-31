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
    Foam::NSRDSfunc14

Description
    NSRDS function number 114

    Source:
    @verbatim
                      NSRDS - AICHE
                 Data Compilation Tables
                    of Properties of
                     Pure Compounds

        Design Institute for Physical Property Data
          American Institute of Chemical Engineers
                  345 East 47th Street
                New York, New York 10017

         National Standard Reference Data System
         American Institute of Chemical Engineers

          T.E. Daubert       -       R.P. Danner

            Department of Chemical Engineering
            The Pennsylvania State University
                University Park, PA 16802
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef NSRDSfunc14_H
#define NSRDSfunc14_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class NSRDSfunc14 Declaration
\*---------------------------------------------------------------------------*/

class NSRDSfunc14
:
    public thermophysicalFunction
{
    // Private data

        // NSRDS function 114 coefficients
        scalar Tc_, a_, b_, c_, d_;


public:

    //- Runtime type information
    TypeName("NSRDSfunc14");


    // Constructors

        //- Construct from components
        NSRDSfunc14(scalar Tc, scalar a, scalar b, scalar c, scalar d)
        :
            Tc_(Tc),
            a_(a),
            b_(b),
            c_(c),
            d_(d)
        {}

        //- Construct from Istream
        NSRDSfunc14(Istream& is)
        :
            Tc_(readScalar(is)),
            a_(readScalar(is)),
            b_(readScalar(is)),
            c_(readScalar(is)),
            d_(readScalar(is))
        {}


    // Member Functions

        //- Evaluate the function and return the result
        scalar f(scalar, scalar T) const
        {
            scalar t = 1.0 - T/Tc_;
            return a_*a_/t + b_ - t*(2.0*a_*c_ + t*(a_*d_ + t*(c_*c_/3.0 + t*(0.5*c_*d_ + 0.2*d_*d_*t))));
        }


        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << Tc_ << token::SPACE
                << a_ << token::SPACE
                << b_ << token::SPACE
                << c_ << token::SPACE
                << d_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const NSRDSfunc14& f)
        {
            f.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
