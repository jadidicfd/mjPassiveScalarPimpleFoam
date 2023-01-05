**************************************MJ NOTE***********************************************************************************

	Version01
	Last modified: March 28, 2021 by MJ
   
	Tested on: OpenFOAM v2012 , v2006
	Note:
        1. This solver MUST be used for Post-processing of turbulance results
         
        2. It is the solver plus an two essential utilities to calculate budgests of the mean turbulent kinetic
            energy equation,turbulent heat flux, and wall heat flux for imcomperssible solver.
            
        3. I have added some comments to the main solver which should make it easy to understand the terms of TKE budgets and heat flux.
         
        4. Both solver and utilities ("mjPostProcessTKEBudgestsLES" and "mjPostProcessMeanWallHeatFluxIncompressible") have been compiled on OF v2012 and v2006.
         
        5. Only some important notes remain:
            5-1-This implementation requires a well established averaged Velocity field (UMean and TMean) to exist in the time folder.
              The process of averaging of TKE budgets happens within the main solver and UMean and TMean are NOT updated within the code.
              So a reasonable averaged field must exist prior to beginning of simulation using "mjPassiveScalarPimpleFoam".
              To provide this averaged field, one can perform simulation for several flow-through time using the native "mjPassiveScalarPimpleFoam" solver.

            5-2-Three of TKE budgets are calculated afterwards the simulation through the "mjPostProcessTKEBudgestsLES" utilitiy provided.
              This utility requires averaged fields, including:

                  -UMean
                  -kMean
                  -UPrime2Mean
                  -turbDiffMean
                  -SGSDiffMean

            If any of the above fields would not exist, the utility will be terminated.

            5-3-turbulent heat flux (UPrime*TPrime) sould be time-averaged

            5-4-wallHeatFlux are calculated afterwards the simulation through the "mjPostProcessMeanWallHeatFluxIncompressible" utilitiy provided.
            This utility requires averaged fields, including:

              -TMean
              -alphatMean

            Moreover, the following properties should be included in transportProperties dictionary:
                    // Laminar Prandtl number
                      dimensionedScalar Pr("Pr", dimless, laminarTransport);
                    // Turbulent Prandtl number
                    dimensionedScalar Prt("Prt", dimless, laminarTransport);
                    // Heat capacity
                    dimensionedScalar Cp0("Cp0", dimSpecificHeatCapacity, laminarTransport);
                    // Fluid density
                    dimensionedScalar rho0("rho0", dimDensity, laminarTransport);

            So make sure to activate averaging process for them in the fieldAverage function object within controlDict.
            It is recommended to add all the defined variables in the main solver to the function object.

            6. If you would like to use the solver in your work, please cite the following publication:
                  -Asgari, E., and Tadjfar, M., 2017, “Assessment of Four Inﬂow Conditions on Large-Eddy Simulation of a Gently Curved Backward-Facing Step,” J. Turbul.,
                  18(1), pp. 61–86.
                  -Mohammad.jadidi@manchestr.uk.ac
				  
****************************************************************************************************************************************
	Version02
	Last modified: April 15, 2021 by MJ
	Tested on: OpenFOAM v2012 , v2006
  
		1. createFiels modified. pMean and TMean-->         
		IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
		
		2. For savaing the UMean, TMean and pMean during the TKE budget production, new varialbles have been defined:
				UMeanMJ
				TMeanMJ
				pMeanMj

****************************************************************************************************************************************
	Version03
	Last modified: April 18, 2021 by MJ
	Tested on: OpenFOAM v2012 , v2006		       

		1. calculate the SGS turbulent kinetic energy 
		2. calculate the total (SGS + resolved) turbulent kinetic energy 
		2. calculate turbulent dissipation rate at runtime during LES simulations. 
		3. It can be easily extended to include all the terms of the turbulent kinetic energy budget.

	    4. The LES resolution index is also calculated at runtime.

				[LES_{ResIndex}=\frac{k_{Res}}{k_{Res}+k_{SGS}}]

	Please note that the utlity assumes that the fieldAverage utility is used to calculate UMean.
	This is then used to calculate the fluctuating velocity vector UPrime within tkeBudget.H as UPrime = U - UMean.
	
****************************************************************************************************************************************

****************************************************************************************************************************************
	Version04
	Last modified: April 18, 2021 by MJ
	Tested on: OpenFOAM v2012 , v2006		       

		1. cleanup of the code
	
****************************************************************************************************************************************

****************************************************************************************************************************************
	Version05
	Last modified: April 18, 2021 by MJ
	Tested on: OpenFOAM v2012 , v2006		       

		1. Intergral length scale based on TKE and dissipation is calculated
		2. Kolmogorov Length Scale is calclulated
		
	
****************************************************************************************************************************************

