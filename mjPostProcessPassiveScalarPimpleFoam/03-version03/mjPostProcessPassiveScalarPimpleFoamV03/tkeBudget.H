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
