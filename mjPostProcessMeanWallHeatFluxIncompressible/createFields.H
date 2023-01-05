    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    Info << "Reading field QMean\n" << endl;
    volScalarField QMean
    (
        IOobject
        (
            "QMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("QMean", dimMass/(dimLength*dimTime*dimTime*dimTime),0)
    );

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );



    Info << "Reading field TMean\n" << endl;
    volScalarField TMean
    (
        IOobject
        (
            "TMean",
            runTime.timeName(),
            mesh,
	          IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );




        Info << "Reading field alphatMean\n" << endl;
        volScalarField alphatMean
        (
          IOobject
            (
                "alphatMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );


    Info << "Creating field alphaEffMean\n" << endl;
    volScalarField alphaEffMean
    (
        IOobject
        (
            "alphaEffMean",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaEffMean", dimLength*dimLength/dimTime,0)
    );



    Info << "Creating field agradTMeanMJ\n" << endl;
    surfaceScalarField gradTMeanMJ
    (
        IOobject
        (
            "gradTMeanMJ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("gradTMeanMJ", dimTemperature/dimLength,scalar(0))
    );


	#include "createPhi.H"

//    label pRefCell = 0;
//    scalar pRefValue = 0.0;
//    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    //Correct boundaries only after all fields are loaded in
    T.correctBoundaryConditions();
    TMean.correctBoundaryConditions();
