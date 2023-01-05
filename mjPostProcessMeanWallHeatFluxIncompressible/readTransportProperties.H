   
    // Laminar Prandtl number
    dimensionedScalar Pr("Pr", dimless, laminarTransport);

    // Turbulent Prandtl number
    dimensionedScalar Prt("Prt", dimless, laminarTransport);
    
    // Heat capacity
    dimensionedScalar Cp0("Cp0", dimSpecificHeatCapacity, laminarTransport);
    
    // Fluid density
    dimensionedScalar rho0("rho0", dimDensity, laminarTransport);
    

