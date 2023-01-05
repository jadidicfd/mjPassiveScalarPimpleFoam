/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    mjPostProcessPassiveScalarPimpleFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

**************************************MJ NOTE***********************************************************************************

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

              5-4- If you would like to use the solver in your work, please cite the following publication:
                  -Asgari, E., and Tadjfar, M., 2017, “Assessment of Four Inﬂow Conditions on Large-Eddy Simulation of a Gently Curved Backward-Facing Step,” J. Turbul.,
                  18(1), pp. 61–86.
                  -Mohammad.jadidi@manchestr.uk.ac


****************************************************************************************************************************************

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }


        #include "TEqn.H"


        const volScalarField nuLam(turbulence->nu());

		    const volSymmTensorField G(turbulence->R());

		    const volScalarField nut(turbulence->nut());



        UPrime = (U-UMean);         //Resolved Velocity Fluctuation Vector

        pPrime = (p-pMean);              //Resolved Pressure Fluctuation

        B = -2.0*nut*symm(fvc::grad(U));

        turbDiff = -0.5*(UPrime*magSqr(UPrime));      //Turbulent Diffusion Term--divergence operator will be applied afterwards

        pressDiff = -UPrime & fvc::grad(pPrime);           //Pressure Diffusion Term

        SGSstrainTensor = symm(fvc::grad(UPrime));       //Tensor of Strain Rate of Resolved Fluctuations

        viscDiss = -2*nuLam*(SGSstrainTensor && SGSstrainTensor); //Viscous Dissipation of Resolved Fluctuations

        SGSDiff = -UPrime & G;                        //SGS Diffusion Term--divergence operator will be applied afterwards

        SGSDiss = B && SGSstrainTensor;                  //SGS Dissipation Term


        TPrime = (T-TMean);                             //Resolved Temperature Fluctuation Scalar
        turbHeatFlux = UPrime*TPrime;                   //Resolved Turbulent Heat FLux Vector


        runTime.write();
        runTime.printExecutionTime(Info);

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
