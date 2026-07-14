	//- Read the comfortFoam configuration dictionary
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

	//- Personal parameters and ISO 7730:2025 preprocessing options

        // Clothing insulation [clo]
        const scalar clo (readScalar(comfortFoamDict.lookup("clo")));

        // Metabolic rate [met]
        const scalar met (readScalar(comfortFoamDict.lookup("met")));

        // External work [met]
        const scalar wme (readScalar(comfortFoamDict.lookup("wme")));

        // Fallback relative humidity [%]
        const scalar RH1 (readScalar(comfortFoamDict.lookup("RH")));

        // Exactly one of these should be true for ISO 7730:2025 results.
        const bool applyDynamicClothing
        (
            comfortFoamDict.lookupOrDefault<bool>("applyDynamicClothing", false)
        );
        const bool clothingAlreadyResultant
        (
            comfortFoamDict.lookupOrDefault<bool>("clothingAlreadyResultant", false)
        );

        // -1 selects the activity-based value when dynamic correction is enabled.
        const scalar configuredWalkingSpeed
        (
            comfortFoamDict.lookupOrDefault<scalar>("walkingSpeed", -1.0)
        );

        // Optional constant override; a field with this name has priority.
        const scalar configuredRelativeAirVelocity
        (
            comfortFoamDict.lookupOrDefault<scalar>("relativeAirVelocity", -1.0)
        );

        // Optional direct water-vapour pressure [Pa]. RH is used when negative.
        const scalar configuredWaterVapourPressure
        (
            comfortFoamDict.lookupOrDefault<scalar>("waterVapourPressure", -1.0)
        );

        // Optional local-discomfort inputs. Negative values mean unavailable.
        const scalar configuredVerticalTemperatureDifference
        (
            comfortFoamDict.lookupOrDefault<scalar>
            (
                "verticalAirTemperatureDifference",
                -1.0
            )
        );
        const scalar configuredFloorTemperature
        (
            comfortFoamDict.lookupOrDefault<scalar>("floorTemperature", -1.0)
        );
        const scalar configuredRadiantTemperatureAsymmetry
        (
            comfortFoamDict.lookupOrDefault<scalar>
            (
                "radiantTemperatureAsymmetry",
                -1.0
            )
        );
        const word radiantAsymmetryType
        (
            comfortFoamDict.lookupOrDefault<word>
            (
                "radiantAsymmetryType",
                "warmCeiling"
            )
        );


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

    IOobject MRTHeader
    (
        "MRT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject relativeAirVelocityHeader
    (
        "relativeAirVelocity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject waterVapourPressureHeader
    (
        "waterVapourPressure",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject verticalTemperatureDifferenceHeader
    (
        "verticalAirTemperatureDifference",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject floorTemperatureHeader
    (
        "floorTemperature",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    IOobject radiantTemperatureAsymmetryHeader
    (
        "radiantTemperatureAsymmetry",
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
        dimensionedScalar("DR", dimensionSet(0,0,0,0,0,0,0), 0.0),
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

    volScalarField PDVertical
    (
        IOobject
        (
            "PDVertical",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PDVertical", dimless, -1.0),
        "zeroGradient"
    );

    volScalarField PDFloor
    (
        IOobject
        (
            "PDFloor",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PDFloor", dimless, -1.0),
        "zeroGradient"
    );

    volScalarField PDRadiantAsymmetry
    (
        IOobject
        (
            "PDRadiantAsymmetry",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PDRadiantAsymmetry", dimless, -1.0),
        "zeroGradient"
    );

    volScalarField ISO7730Valid
    (
        IOobject
        (
            "ISO7730Valid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ISO7730Valid", dimless, -1.0),
        "zeroGradient"
    );

    volScalarField ISO7730WholeBodyCategory
    (
        IOobject
        (
            "ISO7730WholeBodyCategory",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ISO7730WholeBodyCategory", dimless, -1.0),
        "zeroGradient"
    );

    volScalarField ISO7730Category
    (
        IOobject
        (
            "ISO7730Category",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("ISO7730Category", dimless, -1.0),
        "zeroGradient"
    );
