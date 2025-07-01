Info<< nl << "Calculating specific catch ratio fields\n" << endl;
for (int phase_no = 0; phase_no < phases.size(); phase_no++)
{
	string phase_scrname = "scr";
	string phase_no_str;         
	std::stringstream convert;
	convert << phase_no+1;
	phase_no_str = convert.str();
	phase_scrname.append(phase_no_str);
	
	tmp<surfaceScalarField> tnormalvel
	(
		new surfaceScalarField
		(
			IOobject
			(
				"normalvel",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mag( (mesh.Sf()/mesh.magSf()) & fvc::interpolate(Urain[phase_no]) )
		)
	);
	const surfaceScalarField& normalvel = tnormalvel();

	scr.set
	(
		phase_no,
		new volScalarField
		(
			IOobject
			(
				phase_scrname,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			(fvc::average(normalvel) * alpharain[phase_no]) * ((3600*1E3)/(Rh*phases[phase_no][1]))
		)
	);
	scr[phase_no].write();
}

Info<< nl << "Calculating catch ratio field\n" << endl;

volScalarField gcr //initialize global catch ratio field
(
	IOobject
	(
		"gcr",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	scr[0]*phases[0][1]
);

for (int phase_no = 1; phase_no < phases.size(); phase_no++)
{
	gcr += scr[phase_no]*phases[phase_no][1];
}
gcr.write();
