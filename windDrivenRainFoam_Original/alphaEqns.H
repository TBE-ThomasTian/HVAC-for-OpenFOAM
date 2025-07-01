fvScalarMatrix alphaEqn
(
  fvm::ddt(alpharain[phase_no]) 
+ fvm::div(phirain[phase_no],alpharain[phase_no],"div(phirain,alpharain)")  
);

alphaEqn.relax();

alphaEqn.solve();

alpharain[phase_no].correctBoundaryConditions();
