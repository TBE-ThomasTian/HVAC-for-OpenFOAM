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
        
        //- Mass flow of air m3
    //    const scalar sumZuluft (readScalar(comfortFoamDict.lookup("Zuluft")));

		// Room humidity 
        const scalar RH1 (readScalar(comfortFoamDict.lookup("RH")));  


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
        dimensionedScalar("G", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volScalarField w
    (
        IOobject
        (
            "w",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("w", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    volScalarField RH
    (
        IOobject
        (
            "RH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("RH", dimensionSet(0,0,0,0,0,0,0), 0.0)
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
        dimensionedScalar("TOp", dimensionSet(0,0,0,1,0,0,0), 0.0)
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
        dimensionedScalar("DR", dimensionSet(0,1,-1,0,0,0,0), 0.0)
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
		dimensionedScalar("PMV", dimensionSet(0,0,0,0,0,0,0), 0.0)
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
        dimensionedScalar("PPD", dimensionSet(0,0,0,0,0,0,0), 0.0)
	);

	bool RHswitch = false;

