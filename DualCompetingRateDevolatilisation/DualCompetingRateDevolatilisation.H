/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::DualCompetingRateDevolatilisation

Group
    grpLagrangianIntermediateDevolatilisationSubModels

Description
    Single kinetic rate devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = A1.exp(- E/R.T)

    Where:
        kappa = rate constant
        A1    = activation energy (user input)
        E     = pre-exponential factor (user input)
        R     = universal gas constant 
        T     = temperature

    Usage:

        dualCompetingRateDevolatilisationCoeffs
        {
            volatileData
            (
                (CH4     0.l5	1e3     1e5
			 1.00	1e5	1e5)   // (name a1 A1 E1 a2 A2 E2)
                (CO2     0.l5	1e3	1e5
			 1.00	1e5	1e5)   // (name a1 A1 E1 a2 A2 E2)
            );

            volatileResidualCoeff 1e-6;
        }
	
    NOTE:
	
	R uses kmol insteasd of mol. 
	
\*---------------------------------------------------------------------------*/

#ifndef DualCompetingRateDevolatilisation_H
#define DualCompetingRateDevolatilisation_H

#include "DevolatilisationModel.H"

// custom include



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class DualCompetingRateDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DualCompetingRateDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
        // Private data

            //- Specie name
            word name_;

	    //- Stoichiometric coefficient - 1
	    scalar a1_;

            //- Pre-exponential factor - 1
            scalar A1_;

            //- Activiation Energy - 1
            scalar E1_;

	    //- Stoichiometric coefficient - 1
	    scalar a2_;

            //- Pre-exponential factor - 2
            scalar A2_;

            //- Activiation Energy - 2
            scalar E2_;

	    // Addressing

            //- Ash positions in global/local lists
            label ashLocalId_;

    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
		a1_(0.0),
                A1_(0.0),
                E1_(0.0),
		a2_(0.0),
		A2_(0.0),
                E2_(0.0)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
		a1_(readScalar(is)),
                A1_(readScalar(is)),
                E2_(readScalar(is)),
		a2_(readScalar(is)),
                A2_(readScalar(is)),
                E2_(readScalar(is))
            {}

            //- Construct as copy
            volatileData(const volatileData& vd)
            :
                name_(vd.name_),
		a1_(vd.a1_),
                A1_(vd.A1_),
                E1_(vd.E1_),
		a2_(vd.a2_),
                A2_(vd.A2_),
                E2_(vd.E2_)
            {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the alpha - 1
                scalar a1() const
                {
                    return a1_;
                }

                //- Return const access to the pre-exponential factor - 1
                scalar A1() const
                {
                    return A1_;
                }

                //- Return const access to the activation energy - 1
                scalar E1() const
                {
                    return E1_;
                }

                //- Return const access to the alpha - 2
                scalar a2() const
                {
                    return a2_;
                }

                //- Return const access to the pre-exponential factor - 2
                scalar A2() const
                {
                    return A2_;
                }

                //- Return const access to the activation energy - 2
                scalar E2() const
                {
                    return E2_;
                }

        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.a1_ >> vd.A1_ >> vd.E1_ >> vd.a2_ >> vd.A2_ >> vd.E2_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
		    << vd.a1_ << token::SPACE
                    << vd.A1_ << token::SPACE
                    << vd.E1_ << token::SPACE
		    << vd.a2_ << token::SPACE
                    << vd.A2_ << token::SPACE	
                    << vd.E2_ << token::END_LIST;

                return os;
            }
    };


    // Private data

        // Model constants

            //- List of volatile data - (name A1&2 E1&2)
            List<volatileData> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;

	   //- Dynamic Volatile Switch 
//	    const scalar dynamicVolatiles_; 

public:

    //- Runtime type information
    TypeName("dualCompetingRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        DualCompetingRateDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        DualCompetingRateDevolatilisation
        (
            const DualCompetingRateDevolatilisation<CloudType>& dm
        );

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DevolatilisationModel<CloudType>>
            (
                new DualCompetingRateDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DualCompetingRateDevolatilisation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,		// timestep
            const scalar age,		// age
            const scalar mass0,		// mass (initial on injection)
            const scalar mass,		// mass
            const scalar T,		// temperature
            const scalarField& YGasEff,	// gas component mass fractions
            const scalarField& YLiquidEff, // liquid component mass fractions
            const scalarField& YSolidEff, // solid component mass fractions
            label& canCombust,		// 'can combust' flag
            scalarField& dMassDV	// mass transfer - local to particle 
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DualCompetingRateDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
