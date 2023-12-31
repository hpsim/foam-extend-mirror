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

#include "spray.H"
#include "breakupModel.H"
#include "collisionModel.H"
#include "dispersionModel.H"
#include "injectorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void spray::inject()
{
    scalar time = runTime_.value();
    scalar time0 = time0_;

    // Inject the parcels for each injector sequentially
    forAll(injectors_, i)
    {
        autoPtr<injectorType>& it = injectors()[i].properties();
        if (!it->pressureIndependentVelocity())
        {
            scalar referencePressure = p().average().value();
            it->correctProfiles(fuels(), referencePressure);
        }

        const label nHoles = it->nHoles();

        // parcels have the same mass during a timestep
        scalar mass = it->mass(time0, time, twoD_, angleOfWedge_);

        label Np = it->nParcelsToInject(time0, time);

        if (mass > 0)
        {
            Np = max(1, Np);
            scalar mp = mass/Np/nHoles;

            // constT is only larger than zero for the first
            // part of the injection
            scalar constT = max(0.0, it->tsoi() - time0);

            // deltaT is the duration of injection during this timestep
            scalar deltaT = min
            (
                runTime_.deltaT().value(),
                min
                (
                    time - it->tsoi(),
                    it->teoi() - time0
                )
            );

            for(label j=0; j<Np; j++)
            {
                // calculate the time of injection for parcel 'j'
                scalar toi = time0 + constT + deltaT*j/scalar(Np);

                for(label n=0; n<nHoles; n++)
                {

                    // calculate the velocity of the injected parcel
                    vector injectionPosition = it->position
                    (
                        n,
                        toi,
                        twoD_,
                        angleOfWedge_,
                        axisOfSymmetry_,
                        axisOfWedge_,
                        axisOfWedgeNormal_,
                        rndGen_
                    );

                    scalar diameter = injection().d0(i, toi);
                    vector direction = injection().direction(i, n, toi, diameter);
                    vector U = injection().velocity(i, toi)*direction;

                    scalar symComponent = direction & axisOfSymmetry_;
                    vector normal = direction - symComponent*axisOfSymmetry_;
                    normal /= mag(normal);

                    // should be set from dict or model
                    scalar deviation = breakup().y0();
                    scalar ddev = breakup().yDot0();

                    label injectorCell = mesh_.findCell(injectionPosition);

#                   include "findInjectorCell.H"

                    if (injectorCell >= 0)
                    {
                        scalar liquidCore = 1.0;

                        // construct the parcel that is to be injected

                        parcel* pPtr = new parcel
                        (
                            *this,
                            injectionPosition,
                            injectorCell,
                            normal,
                            diameter,
                            it->T(toi),
                            mp,
                            deviation,
                            ddev,
                            0.0,
                            0.0,
                            0.0,
                            liquidCore,
                            scalar(i),
                            U,
                            vector::zero,
                            it->X(),
                            fuels_->components()
                        );

                        injectedLiquidKE_ += 0.5*pPtr->m()*magSqr(U);

                        scalar dt = time - toi;

                        pPtr->stepFraction() =
                            (runTime_.deltaT().value() - dt)
                           /runTime_.deltaT().value();

                        bool keepParcel = pPtr->move(*this);

                        if (keepParcel)
                        {
                            addParticle(pPtr);
                        }
                        else
                        {
                            delete pPtr;
                        }
                    } // if (injectorCell....
                } // for(label n=0...
            } // for(label j=0....
        } // if (mass>0)...
    } // forAll(injectors)...

    time0_ = time;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
