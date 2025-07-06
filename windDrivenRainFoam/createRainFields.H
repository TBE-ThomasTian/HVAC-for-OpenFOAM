PtrList<volVectorField> Urain;
Urain.setSize(phases.size());
PtrList<surfaceScalarField> phirain;
phirain.setSize(phases.size());
PtrList<volScalarField> alpharain;
alpharain.setSize(phases.size());
PtrList<volScalarField> scr;
scr.setSize(phases.size());

std::stringstream converter; 	
for (int phase_no = 0; phase_no < phases.size(); phase_no++)
{
	string phase_Uname = "U";
	string phase_phiname = "phi";
	string phase_alphaname = "alpha";
	
	string phase_no_str;
	std::stringstream convert;
	convert << phase_no+1;
	phase_no_str = convert.str();
	
	phase_Uname.append(phase_no_str);
	phase_phiname.append(phase_no_str);
	phase_alphaname.append(phase_no_str);

	tmp<volVectorField> Utemp(new volVectorField
	(
		IOobject
		(
			phase_Uname,
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	));
	Urain.set(phase_no,Utemp);

	tmp<surfaceScalarField> phitemp(new surfaceScalarField
	(
		IOobject
		(
			phase_phiname,
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		linearInterpolate(Urain[phase_no]) & mesh.Sf()
	));
	phirain.set(phase_no,phitemp);

	tmp<volScalarField> alphatemp(new volScalarField
	(
		IOobject
		(
			phase_alphaname,
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	));
	alpharain.set(phase_no,alphatemp);
	
}

// Create fields for turbulent dispersion (always created, but only used if solveTD)
PtrList<volScalarField> Ctrain;
Ctrain.setSize(phases.size());	

for (int phase_no = 0; phase_no < phases.size(); phase_no++)
{
    string phase_Ctname = "Ct";
    string phase_no_str;         
    std::stringstream convert;
    convert << phase_no+1;
    phase_no_str = convert.str();
    phase_Ctname.append(phase_no_str);
    
    Ctrain.set
    (
        phase_no,
        new volScalarField
        (
            IOobject
            (
                phase_Ctname,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(phase_Ctname, dimless, 0)
        )
    );
}

volScalarField nutrain
(
    IOobject
    (
        "nutrain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("nutrain", dimensionSet(0,2,-1,0,0,0,0), 1)
);

volScalarField tp
(
    IOobject
    (
        "tp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tp", dimTime, 0)
);

volScalarField tfl
(
    IOobject
    (
        "tfl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tfl", dimTime, 0)
);
