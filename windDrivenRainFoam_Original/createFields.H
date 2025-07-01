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
dimensionedScalar Rh
(
	transportProperties.lookup("Rh")
);
dimensionedScalar temp
(
	transportProperties.lookup("temp")
);	
dimensionedScalar rhoa
(
	"rhoa",
	dimensionSet(1,-3,0,0,0,0,0),
	scalar(0) 
);
dimensionedScalar rhop
(
	"rhop",
	dimensionSet(1,-3,0,0,0,0,0),
	scalar(0) 
);
dimensionedScalar mua
(
	"mua",
	dimensionSet(1,-1,-1,0,0,0,0),
	scalar(0) 
);
scalar scalingFactor
(
	transportProperties.lookupOrDefault<scalar>("scalingFactor", 1.0)
);
const bool solveTD = transportProperties.lookupOrDefault("solveTD", false);	
List<List<scalar> > phases
(
transportProperties.lookup("phases")
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
	IOobject
	(
		"U",
		"0",//runTime.timeName(), //If read from the latest timedir, it will be scaled each time the simulation is restarted!
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	),
	mesh
);
U = U*scalingFactor;
U.correctBoundaryConditions();

surfaceScalarField phiwind
(
	IOobject
	(
		"phiwind",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	linearInterpolate(U) & mesh.Sf()
);

volScalarField Re
(
	IOobject
	(
		"Re",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("Re",dimensionSet(0,0,0,0,0,0,0),500)
);

volScalarField CdRe
(
	IOobject
	(
		"CdRe",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("CdRe",dimensionSet(0,0,0,0,0,0,0),250)
);
