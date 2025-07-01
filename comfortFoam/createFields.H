	//- Read a dictionary to configure comfortTool
	IOdictionary comfortFoamDict
	(
    	IOobject
		(
			"comfortFoamDict",
			runTime.time().constant(),
			runTime,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	//- Define values for the comfortTool from dictionary

        // Clothing
        const scalar clo (readScalar(comfortFoamDict.lookup("clo")));

	    // Metabolic rate
        const scalar met (readScalar(comfortFoamDict.lookup("met")));

        // External work
        const scalar wme (readScalar(comfortFoamDict.lookup("wme")));

       	// Room humidity
        const scalar RH1 (readScalar(comfortFoamDict.lookup("RH")));


	IOobject THeader
	(
		"T",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ
	);

	IOobject kHeader
	(
		"k",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
	);

	IOobject epsilonHeader
	(
		"epsilon",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
	);

	IOobject omegaHeader
	(
		"omega",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
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
        dimensionedScalar("G", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
    );

    // Additional radiation fields for different radiation models
    IOobject qrHeader
    (
        "qr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject IDefaultHeader
    (
        "IDefault",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    volScalarField relHum
    (
        IOobject
        (
            "relHum",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("relHum", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
    );

    volScalarField RH
    (
        IOobject
        (
            "RH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("RH", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
    );

    IOobject UHeader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject p_rghHeader
    (
        "p_rgh",
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

    volScalarField DR
    (
        IOobject
        (
            "DR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("DR", dimensionSet(0,1,-1,0,0,0,0), 0.0),
        "zeroGradient"
	);

    volScalarField PMV
    (
        IOobject
        (
            "PMV",
            runTime.timeName(),
            mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
        ),
		mesh,
		dimensionedScalar("PMV", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
	);

	volScalarField PPD
    (
        IOobject
        (
            "PPD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PPD", dimensionSet(0,0,0,0,0,0,0), 0.0),
        "zeroGradient"
	);

