	IOobject THeader
	(
		"T",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ
	);

    volScalarField G
    (
        IOobject
        (
        	"G",
        	runTime.timeName(),
        	mesh,
        	IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    	),
        mesh,
        dimensionedScalar("G", dimensionSet(1,0,-3,0,0,0,0), 0.0),
        "zeroGradient"
    );


    IOobject UHeader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );


    volScalarField TOp
    (
        IOobject
        (
            "TOp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TOp", dimensionSet(0,0,0,1,0,0,0), 0.0),
        "zeroGradient"
    );

	volScalarField ASHRAELevel80
    (
        IOobject
        (
            "ASHRAELevel80",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ASHRAELevel80", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
	);

	volScalarField ASHRAELevel90
    (
        IOobject
        (
            "ASHRAELevel90",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ASHRAELevel90", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
	);

