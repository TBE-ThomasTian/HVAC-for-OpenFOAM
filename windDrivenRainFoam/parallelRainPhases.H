// Parallel processing of rain phases using OpenMP
// This allows multiple rain phases to be computed simultaneously on shared memory systems

#ifdef _OPENMP
    #include <omp.h>
    Info<< "OpenMP parallelization enabled with " << omp_get_max_threads() << " threads" << endl;
#endif

// Function to process a single rain phase
void processRainPhase
(
    const int phase_no,
    const scalar dp,
    const scalar rhoa,
    const scalar rhop,
    const scalar mua,
    const volVectorField& U,
    const volVectorField& g,
    const bool solveTD,
    const autoPtr<incompressible::turbulenceModel>& turbulence,
    volScalarField& Re,
    volScalarField& CdRe,
    volScalarField& tp,
    volScalarField& tfl,
    volScalarField& nutrain,
    PtrList<volVectorField>& Urain,
    PtrList<surfaceScalarField>& phirain,
    PtrList<volScalarField>& alpharain,
    PtrList<volScalarField>& Ctrain,
    const fvMesh& mesh
)
{
    // Calculate relative velocity
    volVectorField Ur = U - Urain[phase_no];
    volScalarField magUr = mag(Ur);
    
    // Reynolds number and drag coefficient
    Re = (magUr*dp*rhoa)/mua;
    CdRe = GET_CdRe(Re);
    CdRe.correctBoundaryConditions();
    
    // Turbulent dispersion if enabled
    if (solveTD && turbulence.valid())
    { 
        // Calculate turbulent time scales
        tmp<volScalarField> tk = turbulence->k();
        tmp<volScalarField> tepsilon = turbulence->epsilon();
        tfl = 0.2*(tk()/tepsilon());
        tp = (4*rhop*dp*dp)/(3*mua*CdRe);
        Ctrain[phase_no] = sqrt( tfl/(tfl+tp) );
        nutrain = turbulence->nut()()*sqr(Ctrain[phase_no]);
    }
    else
    {
        nutrain = dimensionedScalar("zero", nutrain.dimensions(), 0.0);
    }
    
    // Include momentum equation
    #include "UEqnSinglePhase.H"
}

// Parallel phase processing structure
struct PhaseComputeData
{
    int phase_no;
    scalar dp;
    volScalarField Re;
    volScalarField CdRe;
    volScalarField tp;
    volScalarField tfl;
    volScalarField nutrain;
    
    PhaseComputeData(const fvMesh& mesh, int pn, scalar diameter)
    :
        phase_no(pn),
        dp(diameter),
        Re
        (
            IOobject
            (
                "Re_" + std::to_string(pn),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Re", dimless, 500)
        ),
        CdRe
        (
            IOobject
            (
                "CdRe_" + std::to_string(pn),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("CdRe", dimless, 250)
        ),
        tp
        (
            IOobject
            (
                "tp_" + std::to_string(pn),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("tp", dimTime, 0)
        ),
        tfl
        (
            IOobject
            (
                "tfl_" + std::to_string(pn),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("tfl", dimTime, 0)
        ),
        nutrain
        (
            IOobject
            (
                "nutrain_" + std::to_string(pn),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("nutrain", dimensionSet(0,2,-1,0,0,0,0), 1)
        )
    {}
};