//part1

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

//part2

    volScalarField kSGS = turbulence->k(); //Instantaneous subgrid-scale tke

    volScalarField kRes = 0.5*(UPrime & UPrime); //Instantaneous resolved tke

    kTot = kRes + kSGS;

    volScalarField epsilonSGS = turbulence->epsilon();

    volSymmTensorField SRes = symm(fvc::grad(UPrime)); //Instantaneous resolved fluctuating rate of strain
    volScalarField epsilonRes=-2*nuLam*(SRes && SRes); //Instantaneous resolved dissipation rate

    epsilonTot = mag(epsilonRes + epsilonSGS);

    LESResIndex = kRes/max(kSmall,kTot); //Instantaneous ratio of resolved to total (Resolved + SGS) tke

                  Info << "ESTIMATED LESResIndex min, max, avg: " << min(LESResIndex).value()
              					    << ", " << max(LESResIndex).value()
              					    << ", " << LESResIndex.weightedAverage(mesh.V()).value()
              					    << endl;

//part3

    volScalarField integralLengthScale=Foam::pow(kTot, 3.0 / 2.0) / epsilonTot;

//part4

    TPrime = (T-TMean);                             //Resolved Temperature Fluctuation Scalar
    turbHeatFlux = UPrime*TPrime;                   //Resolved Turbulent Heat FLux Vector
