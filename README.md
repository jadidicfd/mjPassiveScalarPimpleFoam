# mjPassiveScalarPimpleFoam
Transient solver for incompressible turbulent flow with a passive scalar
Description: Transient solver for incompressible, turbulent flow of Newtonian fluids. Heading Solver details: The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm. Sub-models include: - turbulence modelling, i.e. laminar, RAS or LES; - run-time selectable MRF and finite volume options, e.g. explicit porosity. Solver capability: 1. This solver can also be used for Post-processing of turbulence results. 2. It is the solver plus two essential utilities to calculate budgets of the mean turbulent kinetic energy equation, turbulent heat flux, and wall heat flux for the incompressible solver. 3. Comments have been added to the main solver which should make it easy to understand the terms of TKE budgets and heat flux. 4. Both solver and utilities ("mjPostProcessTKEBudgestsLES" and "mjPostProcessMeanWallHeatFluxIncompressible") have been compiled on OF v2012 and v2006. 5-1- The implementation requires a well-established averaged Velocity field (UMean and TMean) to exist in the time folder. The process of averaging of TKE budgets happens within the main solver and UMean and TMean are NOT updated within the code. So a reasonable averaged field must exist prior to beginning of simulation using "mjPassiveScalarPimpleFoam". To provide this averaged field, one can perform simulation for several flow-through time using the native "mjPassiveScalarPimpleFoam" solver. 5-2-Three of TKE budgets are calculated afterward the simulation through the "mjPostProcessTKEBudgestsLES" utility provided. This utility requires averaged fields, including: -UMean -kMean -UPrime2Mean -turbDiffMean -SGSDiffMean If any of the above fields would not exist, the utility will be terminated. 5-3-turbulent heat flux (UPrime*TPrime) should be time-averaged. 5-4-wallHeatFlux are calculated afterward the simulation through the "mjPostProcessMeanWallHeatFluxIncompressible" utility provided. This utility requires averaged fields, including: -TMean -alphatMean NOTE: The solver is now under validation procedure.
