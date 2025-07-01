// Adaptive time stepping for rain phases

scalar maxAlphaCo = 0.0;
scalar maxUrainCo = 0.0;

if (adjustTimeStep)
{
    // Calculate max Courant numbers for all rain phases
    for (int phase_no = 0; phase_no < phases.size(); phase_no++)
    {
        // Volume fraction Courant number
        scalar alphaCoNum = 0.0;
        if (alpharain[phase_no].dimensions() == dimless)
        {
            surfaceScalarField phiAlpha = fvc::interpolate(alpharain[phase_no])*phirain[phase_no];
            alphaCoNum = 0.5*gMax
            (
                fvc::surfaceSum(mag(phiAlpha))().primitiveField()
               /mesh.V().field()
            )*runTime.deltaTValue();
        }
        
        // Velocity-based Courant number
        scalar UrainCoNum = 0.5*gMax
        (
            fvc::surfaceSum(mag(phirain[phase_no]))().primitiveField()
           /mesh.V().field()
        )*runTime.deltaTValue();
        
        maxAlphaCo = max(maxAlphaCo, alphaCoNum);
        maxUrainCo = max(maxUrainCo, UrainCoNum);
        
        Info<< "Phase " << phase_no+1 << ": "
            << "alphaCo = " << alphaCoNum 
            << ", UrainCo = " << UrainCoNum << endl;
    }
    
    Info<< "Max alphaCo = " << maxAlphaCo 
        << ", Max UrainCo = " << maxUrainCo << endl;
    
    // Adjust time step based on maximum Courant numbers
    scalar maxCo = max(maxAlphaCo, maxUrainCo);
    scalar maxDeltaTFact = maxCo > SMALL ? maxCoRain/maxCo : 2.0;
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
    
    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );
    
    Info<< "deltaT = " << runTime.deltaTValue() << endl;
}
