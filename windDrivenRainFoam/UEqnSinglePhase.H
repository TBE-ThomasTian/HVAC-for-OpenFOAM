// Momentum equation for a single rain phase
// This is extracted from UEqns.H to be used in parallel processing

{
    volTensorField Rca(-nutrain*(T(fvc::grad(Urain[phase_no]))));
    if (solveTD && turbulence.valid())
    {
        tmp<volScalarField> tk = turbulence->k();
        Rca = Rca + (2.0/3.0)*sqr(Ctrain[phase_no])*I*tk() - (2.0/3.0)*I*tr(Rca);
    }
    Rca.correctBoundaryConditions();

    surfaceScalarField phiRa
    (
        - fvc::interpolate(nutrain)
        *mesh.magSf()*fvc::snGrad(alpharain[phase_no])/fvc::interpolate(alpharain[phase_no] + scalar(0.001))
    );

    fvVectorMatrix UrainEqn
    (
        fvm::ddt(Urain[phase_no])
        + fvm::div(phirain[phase_no], Urain[phase_no],"div(phirain,Urain)")
        - fvm::Sp(fvc::div(phirain[phase_no]),Urain[phase_no])

        - fvm::laplacian(nutrain, Urain[phase_no])
        + fvc::div(Rca)

        + fvm::div(phiRa, Urain[phase_no],"div(phirain,Urain)")
        - fvm::Sp(fvc::div(phiRa), Urain[phase_no])
        + (fvc::grad(alpharain[phase_no])/(fvc::average(alpharain[phase_no]) + scalar(0.001)) & Rca)
        ==
        g
        + ((3* mua * CdRe)/(4 * rhop * dp * dp))*U
        - fvm::Sp(((3* mua * CdRe)/(4 * rhop * dp * dp)),Urain[phase_no])
    );
    UrainEqn.relax();

    UrainEqn.solve();

    Urain[phase_no].correctBoundaryConditions();

    phirain[phase_no] = fvc::interpolate(Urain[phase_no]) & mesh.Sf();
}