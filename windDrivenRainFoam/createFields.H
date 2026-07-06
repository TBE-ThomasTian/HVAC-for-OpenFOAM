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
dimensionedScalar Rh(transportProperties.lookup("Rh"));
if (Rh.dimensions() == dimless)
{
	Info<< "Interpreting dimensionless Rh as mm/h and converting to m/s" << endl;
	Rh = dimensionedScalar("Rh", dimVelocity, Rh.value()/(3600.0*1000.0));
}
else if (Rh.dimensions() != dimVelocity)
{
	FatalIOErrorInFunction(transportProperties)
		<< "Rh must be either dimensionless in mm/h or have velocity dimensions [0 1 -1 0 0 0 0]"
		<< exit(FatalIOError);
}

word tempName("temp");
if (!transportProperties.found(tempName) && transportProperties.found("temperature"))
{
	tempName = "temperature";
}
dimensionedScalar temp(transportProperties.lookup(tempName));
if (temp.dimensions() == dimless)
{
	temp = dimensionedScalar("temp", dimTemperature, temp.value());
}
else if (temp.dimensions() != dimTemperature)
{
	FatalIOErrorInFunction(transportProperties)
		<< tempName << " must be either dimensionless in K or have temperature dimensions"
		<< exit(FatalIOError);
}
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
	transportProperties.lookupOrDefault<scalar>
	(
		"scalingFactor",
		transportProperties.lookupOrDefault<scalar>
		(
			"Ustar",
			transportProperties.lookupOrDefault<scalar>("windScaling", 1.0)
		)
	)
);
const bool solveTD = transportProperties.lookupOrDefault
(
	"solveTD",
	transportProperties.lookupOrDefault("turbulentDispersion", false)
);	
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
