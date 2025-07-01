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

	
