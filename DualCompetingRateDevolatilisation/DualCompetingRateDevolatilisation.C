/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DualCompetingRateDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DualCompetingRateDevolatilisation<CloudType>::
DualCompetingRateDevolatilisation
(
    const dictionary& dict,
    CloudType& owner
)
:
    DevolatilisationModel<CloudType>(dict, owner, typeName),
    volatileData_(this->coeffDict().lookup("volatileData")),
    YVolatile0_(volatileData_.size()),
    volatileToGasMap_(volatileData_.size()),
    residualCoeff_(this->coeffDict().getScalar("residualCoeff"))
{
    if (volatileData_.empty())
    {
        WarningInFunction
            << "Devolatilisation model selected, but no volatiles defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating volatile species:" << endl;

        // Determine mapping between active volatiles and cloud gas components
        const label idGas = owner.composition().idGas();
        const scalar YGasTot = owner.composition().YMixture0()[idGas];
        const scalarField& YGas = owner.composition().Y0(idGas);
        forAll(volatileData_, i)
        {
            const word& specieName = volatileData_[i].name();
            const label id = owner.composition().localId(idGas, specieName);
            volatileToGasMap_[i] = id;
            YVolatile0_[i] = YGasTot*YGas[id];

            Info<< "    " << specieName << ": particle mass fraction = "
                << YVolatile0_[i] << endl;
        }
    }
}


template<class CloudType>
Foam::DualCompetingRateDevolatilisation<CloudType>::
DualCompetingRateDevolatilisation
(
    const DualCompetingRateDevolatilisation<CloudType>& dm
)
:
    DevolatilisationModel<CloudType>(dm),
    volatileData_(dm.volatileData_),
    YVolatile0_(dm.YVolatile0_),
    volatileToGasMap_(dm.volatileToGasMap_),
    residualCoeff_(dm.residualCoeff_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DualCompetingRateDevolatilisation<CloudType>::
~DualCompetingRateDevolatilisation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DualCompetingRateDevolatilisation<CloudType>::calculate
(
    const scalar dt,
    const scalar age,
    const scalar mass0,
    const scalar mass,
    const scalar T,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV
) const
{
    bool done = true;
    forAll(volatileData_, i)
    {
        const label id = volatileToGasMap_[i];
        const scalar massVolatile0 = mass0*YVolatile0_[i];
        const scalar massVolatile = mass*YGasEff[id];
		
        //Ash must be declared as second variable in the thermophysical file 
        // Combustion allowed once all volatile components evolved
        done = done && (massVolatile <= residualCoeff_*massVolatile0);

        // Model coefficients
        const scalar a1 = volatileData_[i].a1();
        const scalar A1 = volatileData_[i].A1();
        const scalar E1 = volatileData_[i].E1();
        
        const scalar a2 = volatileData_[i].a2();
        const scalar A2 = volatileData_[i].A2();
        const scalar E2 = volatileData_[i].E2();
        
        // Kinetic Rates
        const scalar k1 = A1*exp(-E1/(RR*T));
        const scalar k2 = A2*exp(-E2/(RR*T)); 
        
        // Mass transferred from particle to carrier gas phase
        dMassDV[id] = min((a1*k1 + a2*k2)*massVolatile*dt, massVolatile);
	}
    if (done && canCombust != -1)
    {
        canCombust = 1;
    }
}


// ************************************************************************* //
