
// The log_e of 10 ¬¬
//float ln10 = 2.302585093; 

// Happel model for filtration
float As_Coeffs[3] = { 8.99202433, -7.40461178, 0.26224388 };

float ppi = 3.14159264; 

// Limit lambda to positive values
dimensionedScalar Lambda0
(
    "Lambda0",
    dimensionSet(0, -1, 0, 0, 0, 0 ,0),
    0.0
);

Info<< "Reading field: Hydraulic head (h) \n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

fvScalarMatrix InitFlowEquation
(
    fvm::laplacian(h)
);
    InitFlowEquation.solve();

Info<< "Reading field: Porosity (theta) \n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField theta0
(
    theta
);

Info<< "Reading field: Suspended clay (Cw) \n" << endl;
volScalarField Cw
(
    IOobject
    (
        "Cw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field: Stopped clay (Cs) \n" << endl;
volScalarField Cs
(
    IOobject
    (
        "Cs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading bed grain size (ds)\n" << endl;
dimensionedScalar ds
(
    transportProperties.lookup("ds")
);

Info<< "Reading clay size (dp)\n" << endl;
dimensionedScalar dp
(
    transportProperties.lookup("dp")
);

Info<< "Reading water density (rho)\n" << endl;
dimensionedScalar rho
(
    transportProperties.lookup("rho")
);

Info<< "Reading clay density (rho_clay)\n" << endl;
dimensionedScalar rho_clay
(
    transportProperties.lookup("rho_clay")
);

Info<< "Reading sand density (rho_sand)\n" << endl;
dimensionedScalar rho_sand
(
    transportProperties.lookup("rho_sand")
);

Info<< "Reading gravity (g)\n" << endl;
dimensionedScalar g
(
    transportProperties.lookup("g")
);

Info<< "Reading dynamic viscosity (mu)\n" << endl;
dimensionedScalar mu
(
    transportProperties.lookup("mu")
);

Info<< "Reading attachement efficiency (alpha)\n" << endl;
dimensionedScalar alpha
(
    transportProperties.lookup("alpha")
);

Info<< "Reading min porosity (thetaMin)\n" << endl;
dimensionedScalar thetaMin
(
    transportProperties.lookup("thetaMin")
);

Info<< "Reading reference porosity for clay fraction (refPoro)\n" << endl;
dimensionedScalar refPoro
(
    transportProperties.lookup("refPoro")
);

Info<< "Reading min Hydraulic Conductivity (Kmin)\n" << endl;
dimensionedScalar Kmin
(
    transportProperties.lookup("Kmin")
);

Info<< "Reading Hamaker constant (Hamaker)\n" << endl;
dimensionedScalar Hamaker
(
    transportProperties.lookup("Hamaker")
);

Info<< "Reading Boltzmann constan (kBoltz)\n" << endl;
dimensionedScalar kBoltz
(
    transportProperties.lookup("kBoltz")
);

Info<< "Reading Temperature (Temp)\n" << endl;
dimensionedScalar Temp
(
    transportProperties.lookup("Temp")
);

Info<< "Read diffusion coefficient (if negative it's calculated)\n" << endl;
dimensionedScalar DUser
(
    transportProperties.lookup("DUser")
);

// Limit lambda to positive values
dimensionedScalar DT0
(
    "DT0",
    dimensionSet(0, 2, -1, 0, 0, 0 ,0),
    0.0
);

dimensionedScalar DT
(
    max(DT0,DUser)
);

Info<< "Read dispersion alpha (zero if no dispersion)\n" << endl;
dimensionedScalar AlphaDisp
(
    transportProperties.lookup("AlphaDisp")
);

Info<< "Read dispersion power \n" << endl;
dimensionedScalar PowerDisp
(
    transportProperties.lookup("PowerDisp")
);

//Dimensionless numbers for filtration model
// Aspect Ratio
dimensionedScalar N_Ratio
(
    dp/ds 
);
Info<< "Aspect Ratio:"<< N_Ratio.value() <<"\n" << endl;

// van der Waals 
dimensionedScalar N_vdWaals 
( 
    Hamaker/(kBoltz*Temp) 
);
Info<< "van der Waals:"<< N_vdWaals.value() <<"\n" << endl;

// Gravitational 
dimensionedScalar N_gravit 
( 
    4.0*ppi/3.0 * 
    pow(dp/2,4) *
    (rho_clay - rho) *
    g / (kBoltz*Temp) 

);
Info<< "Gravitational:"<< N_gravit.value() <<"\n" << endl;

// Underflow (i.e., flow due to channel slope)
Info<< "Reading Slope\n" << endl;
dimensionedVector Slope 
( 
    transportProperties.lookup("Slope")
);


volScalarField percentClay
(
    IOobject
    (
        "percentClay",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


/*
Info<< "log10-Conductivity log(K)\n" << endl;
volScalarField logK
(
    IOobject
    (
        "logK",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


*/


/*Info<< "Make a velocity field U\n" << endl;
volVectorField U
(
    IOobject
    (
    	"U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/


#include "createFvOptions.H"
