/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2412
    \\  /    A nd           | Web:      www.OpenFOAM.org
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "buildingElementFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildingElementFvPatchScalarField::buildingElementFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    useMultiLayer_(false),
    uValue_(0.24),
    layers_(),
    Tlayers_(),
    gValue_(0.0),
    Rsi_(0.13),
    Rse_(0.04),
    tempExt_(293.0),
    tempSky_(283.0),
    solarExt_(0.0),
    hExt_(25.0),
    emissivity_(0.9),
    viewSky_(0.5),
    radField_("none"),
    debugTemp_(false),
    debugFlux_(false),
    debugFace_(-1),
    debugInterval_(60.0),
    lastDebugTime_(0.0),
    qConvInt_(p.size(), 0.0),
    qRadInt_(p.size(), 0.0),
    qCond_(p.size(), 0.0),
    qConvExt_(p.size(), 0.0),
    qRadExt_(p.size(), 0.0),
    qSolar_(p.size(), 0.0),
    qTotal_(p.size(), 0.0)
{
    refValue() = tempExt_;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::buildingElementFvPatchScalarField::buildingElementFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    useMultiLayer_(dict.found("layers")),
    uValue_(0.0),
    layers_(),
    Tlayers_(),
    gValue_(dict.lookupOrDefault<scalar>("gValue", 0.0)),
    Rsi_(dict.lookupOrDefault<scalar>("Rsi", 0.13)),
    Rse_(dict.lookupOrDefault<scalar>("Rse", 0.04)),
    tempExt_(readScalar(dict.lookup("tempExt"))),
    tempSky_(dict.lookupOrDefault<scalar>("tempSky", tempExt_ - 10.0)),
    solarExt_(dict.lookupOrDefault<scalar>("solarExt", 0.0)),
    hExt_(dict.lookupOrDefault<scalar>("hExt", 25.0)),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0.9)),
    viewSky_(dict.lookupOrDefault<scalar>("viewSky", 0.5)),
    radField_(dict.lookupOrDefault<word>("radField", "none")),
    debugTemp_(dict.lookupOrDefault<bool>("debugTemp", false)),
    debugFlux_(dict.lookupOrDefault<bool>("debugFlux", false)),
    debugFace_(dict.lookupOrDefault<label>("debugFace", -1)),
    debugInterval_(dict.lookupOrDefault<scalar>("debugInterval", 60.0)),
    lastDebugTime_(0.0),
    qConvInt_(p.size(), 0.0),
    qRadInt_(p.size(), 0.0),
    qCond_(p.size(), 0.0),
    qConvExt_(p.size(), 0.0),
    qRadExt_(p.size(), 0.0),
    qSolar_(p.size(), 0.0),
    qTotal_(p.size(), 0.0)
{
    if (useMultiLayer_)
    {
        // Read layer properties
        const List<dictionary> layerDicts(dict.lookup("layers"));
        layers_.setSize(layerDicts.size());
        
        forAll(layerDicts, i)
        {
            layers_[i].lambda = readScalar(layerDicts[i].lookup("lambda"));
            layers_[i].thickness = readScalar(layerDicts[i].lookup("thickness"));
            layers_[i].rho = readScalar(layerDicts[i].lookup("rho"));
            layers_[i].cp = readScalar(layerDicts[i].lookup("cp"));
            layers_[i].nNodes = layerDicts[i].lookupOrDefault<label>("nNodes", 5);
        }
        
        // Initialize temperature fields for each layer
        Tlayers_.setSize(layers_.size());
        forAll(layers_, layerI)
        {
            Tlayers_.set
            (
                layerI,
                new scalarField(layers_[layerI].nNodes * p.size(), tempExt_)
            );
        }
        
        // Calculate effective U-value for information
        uValue_ = calculateUValue();
        
        if (debugTemp_ || debugFlux_)
        {
            Info<< "BuildingElement: Multi-layer wall with " << layers_.size() 
                << " layers" << nl
                << "Effective U-value: " << uValue_ << " W/(m²K)" << endl;
        }
    }
    else
    {
        // Simple mode - just read U-value
        uValue_ = readScalar(dict.lookup("uValue"));
    }

    fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }
}


Foam::buildingElementFvPatchScalarField::buildingElementFvPatchScalarField
(
    const buildingElementFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    useMultiLayer_(ptf.useMultiLayer_),
    uValue_(ptf.uValue_),
    layers_(ptf.layers_),
    Tlayers_(),
    gValue_(ptf.gValue_),
    Rsi_(ptf.Rsi_),
    Rse_(ptf.Rse_),
    tempExt_(ptf.tempExt_),
    tempSky_(ptf.tempSky_),
    solarExt_(ptf.solarExt_),
    hExt_(ptf.hExt_),
    emissivity_(ptf.emissivity_),
    viewSky_(ptf.viewSky_),
    radField_(ptf.radField_),
    debugTemp_(ptf.debugTemp_),
    debugFlux_(ptf.debugFlux_),
    debugFace_(ptf.debugFace_),
    debugInterval_(ptf.debugInterval_),
    lastDebugTime_(ptf.lastDebugTime_),
    qConvInt_(ptf.qConvInt_, mapper),
    qRadInt_(ptf.qRadInt_, mapper),
    qCond_(ptf.qCond_, mapper),
    qConvExt_(ptf.qConvExt_, mapper),
    qRadExt_(ptf.qRadExt_, mapper),
    qSolar_(ptf.qSolar_, mapper),
    qTotal_(ptf.qTotal_, mapper)
{
    // Map layer temperatures if using multi-layer
    if (useMultiLayer_)
    {
        Tlayers_.setSize(ptf.Tlayers_.size());
        forAll(Tlayers_, layerI)
        {
            Tlayers_.set
            (
                layerI,
                new scalarField(ptf.Tlayers_[layerI], mapper)
            );
        }
    }
}


Foam::buildingElementFvPatchScalarField::buildingElementFvPatchScalarField
(
    const buildingElementFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    useMultiLayer_(ptf.useMultiLayer_),
    uValue_(ptf.uValue_),
    layers_(ptf.layers_),
    Tlayers_(),
    gValue_(ptf.gValue_),
    Rsi_(ptf.Rsi_),
    Rse_(ptf.Rse_),
    tempExt_(ptf.tempExt_),
    tempSky_(ptf.tempSky_),
    solarExt_(ptf.solarExt_),
    hExt_(ptf.hExt_),
    emissivity_(ptf.emissivity_),
    viewSky_(ptf.viewSky_),
    radField_(ptf.radField_),
    debugTemp_(ptf.debugTemp_),
    debugFlux_(ptf.debugFlux_),
    debugFace_(ptf.debugFace_),
    debugInterval_(ptf.debugInterval_),
    lastDebugTime_(ptf.lastDebugTime_),
    qConvInt_(ptf.qConvInt_),
    qRadInt_(ptf.qRadInt_),
    qCond_(ptf.qCond_),
    qConvExt_(ptf.qConvExt_),
    qRadExt_(ptf.qRadExt_),
    qSolar_(ptf.qSolar_),
    qTotal_(ptf.qTotal_)
{
    // Deep copy layer temperatures
    if (useMultiLayer_)
    {
        Tlayers_.setSize(ptf.Tlayers_.size());
        forAll(Tlayers_, layerI)
        {
            Tlayers_.set
            (
                layerI,
                new scalarField(ptf.Tlayers_[layerI])
            );
        }
    }
}


Foam::buildingElementFvPatchScalarField::buildingElementFvPatchScalarField
(
    const buildingElementFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    useMultiLayer_(ptf.useMultiLayer_),
    uValue_(ptf.uValue_),
    layers_(ptf.layers_),
    Tlayers_(),
    gValue_(ptf.gValue_),
    Rsi_(ptf.Rsi_),
    Rse_(ptf.Rse_),
    tempExt_(ptf.tempExt_),
    tempSky_(ptf.tempSky_),
    solarExt_(ptf.solarExt_),
    hExt_(ptf.hExt_),
    emissivity_(ptf.emissivity_),
    viewSky_(ptf.viewSky_),
    radField_(ptf.radField_),
    debugTemp_(ptf.debugTemp_),
    debugFlux_(ptf.debugFlux_),
    debugFace_(ptf.debugFace_),
    debugInterval_(ptf.debugInterval_),
    lastDebugTime_(ptf.lastDebugTime_),
    qConvInt_(ptf.qConvInt_),
    qRadInt_(ptf.qRadInt_),
    qCond_(ptf.qCond_),
    qConvExt_(ptf.qConvExt_),
    qRadExt_(ptf.qRadExt_),
    qSolar_(ptf.qSolar_),
    qTotal_(ptf.qTotal_)
{
    // Deep copy layer temperatures
    if (useMultiLayer_)
    {
        Tlayers_.setSize(ptf.Tlayers_.size());
        forAll(Tlayers_, layerI)
        {
            Tlayers_.set
            (
                layerI,
                new scalarField(ptf.Tlayers_[layerI])
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::buildingElementFvPatchScalarField::calculateUValue() const
{
    scalar Rtotal = Rsi_;  // Internal surface resistance
    
    forAll(layers_, i)
    {
        Rtotal += layers_[i].thickness / layers_[i].lambda;
    }
    
    Rtotal += Rse_;  // External surface resistance
    
    return 1.0/Rtotal;
}


void Foam::buildingElementFvPatchScalarField::updateLayerTemperatures()
{
    if (!useMultiLayer_) return;
    
    const scalarField Ti = this->patchInternalField();
    const scalar dt = db().time().deltaTValue();
    
    // For each face on the patch
    forAll(Ti, faceI)
    {
        // Update temperatures through all layers
        
        forAll(layers_, layerI)
        {
            const layer& l = layers_[layerI];
            scalarField& Tlayer = Tlayers_[layerI];
            
            const scalar dx = l.thickness / (l.nNodes - 1);
            const scalar alpha = l.lambda / (l.rho * l.cp);
            const scalar Fo = alpha * dt / (dx * dx);  // Fourier number
            
            // Check stability
            if (Fo > 0.5)
            {
                WarningInFunction
                    << "Fourier number " << Fo << " > 0.5 in layer " << layerI
                    << ". Consider reducing time step or increasing nodes." << endl;
            }
            
            // Update interior nodes (explicit scheme)
            for (label n = 1; n < l.nNodes - 1; n++)
            {
                label idx = faceI * l.nNodes + n;
                Tlayer[idx] = Tlayer[idx] + Fo * 
                    (Tlayer[idx+1] - 2*Tlayer[idx] + Tlayer[idx-1]);
            }
            
            // Boundary conditions between layers
            if (layerI == 0)
            {
                // Inner boundary - coupled to room temperature
                label idx = faceI * l.nNodes;
                Tlayer[idx] = Ti[faceI];  // Or use convection BC
            }
            
            if (layerI == layers_.size() - 1)
            {
                // Outer boundary - coupled to external conditions
                label idx = faceI * l.nNodes + (l.nNodes - 1);
                // Apply external BC including radiation to sky
                // This would need iteration to find surface temperature
                Tlayer[idx] = tempExt_;  // Simplified
            }
        }
    }
}


void Foam::buildingElementFvPatchScalarField::calculateHeatFluxes() const
{
    const scalar sigma = physicoChemical::sigma.value();
    
    // Get temperatures
    const scalarField Ti = this->patchInternalField();
    const scalarField& Tw = *this;  // Wall surface temperature
    
    // Get thermal diffusivity and calculate heat flux from OpenFOAM
    const scalarField& deltaCoeffs = patch().deltaCoeffs();
    scalarField kappa(patch().size());
    
    if (db().foundObject<volScalarField>("alphaEff"))
    {
        const fvPatchScalarField& alphaEff = 
            patch().lookupPatchField<volScalarField, scalar>("alphaEff");
        kappa = alphaEff;
    }
    else if (db().foundObject<volScalarField>("kappaEff"))
    {
        const fvPatchScalarField& kappaEff = 
            patch().lookupPatchField<volScalarField, scalar>("kappaEff");
        kappa = kappaEff;
    }
    else
    {
        // Default value for air
        kappa = 0.026;
    }
    
    // Get internal radiation if available
    if (radField_ != "none" && db().foundObject<volScalarField>(radField_))
    {
        const volScalarField& qr = db().lookupObject<volScalarField>(radField_);
        qRadInt_ = qr.boundaryField()[patch().index()];
    }
    else
    {
        qRadInt_ = 0.0;
    }
    
    // Calculate heat fluxes for each face
    forAll(Ti, faceI)
    {
        // OpenFOAM calculates the heat flux from cell to wall as:
        // q = -k * dT/dn ≈ -k * (Tw - Ti) / delta
        // This includes convection through the boundary layer
        qConvInt_[faceI] = kappa[faceI] * deltaCoeffs[faceI] * (Ti[faceI] - Tw[faceI]);
        
        // Calculate wall properties
        scalar R_total = 1.0 / uValue_;
        scalar R_wall = R_total - Rsi_ - Rse_;  // Pure wall resistance
        
        // Heat flux from interior (calculated by OpenFOAM) - removed, not needed with new approach
        
        // For steady state, we need to solve for T_wall_ext such that:
        // q_in = q_through_wall = q_out
        // where q_out = h_ext*(T_wall_ext - T_ext) + radiation
        
        // First approximation: ignore radiation for T_wall_ext calculation
        // q_in = (T_wall_int - T_wall_ext)/R_wall = h_ext*(T_wall_ext - T_ext)
        // Solving for T_wall_ext:
        scalar h_eff = 1.0/R_wall + hExt_;
        scalar Twall_ext = (Tw[faceI]/R_wall + hExt_*tempExt_) / h_eff;
        
        // Conduction through wall
        qCond_[faceI] = (Tw[faceI] - Twall_ext) / R_wall;
        
        // Exterior convection: positive when outside air is warmer than wall
        qConvExt_[faceI] = hExt_ * (tempExt_ - Twall_ext);
        
        // Exterior radiation: NET flux, positive when wall receives more than it emits
        scalar Trad_eff = viewSky_ * tempSky_ + (1 - viewSky_) * tempExt_;
        qRadExt_[faceI] = emissivity_ * sigma * (pow4(Trad_eff) - pow4(Twall_ext));
        
        // Solar gains: always positive (heat into wall)
        qSolar_[faceI] = gValue_ * solarExt_;
        
        // Energy balance check at INNER surface (should be ~0 in steady state)
        // IN: q_conv_int + q_rad_int (from room)
        // OUT: q_through_wall (to outside)
        // 
        // The heat flux through the wall that must be balanced is what the BC actually implements:
        // From updateCoeffs(), the BC enforces: -k*dT/dn = h_eff*(T_wall - T_ref)
        // Where T_ref = (U*T_ext + h_rad*T_rad + q_solar + q_rad_internal) / h_eff
        // 
        // So the actual heat flux from the room is:
        // q_BC = h_eff*(T_wall - T_ref)
        //      = h_eff*T_wall - (U*T_ext + h_rad*T_rad + q_solar + q_rad_internal)
        //      = h_eff*T_wall - U*T_ext - h_rad*T_rad - q_solar - q_rad_internal
        //      = (U + h_rad)*T_wall - U*T_ext - h_rad*T_rad - q_solar - q_rad_internal
        //      = U*(T_wall - T_ext) + h_rad*(T_wall - T_rad) - q_solar - q_rad_internal
        
        // For validation, compute what the BC is actually doing:
        scalar Trad_eff_local = viewSky_ * tempSky_ + (1 - viewSky_) * tempExt_;
        scalar Tmean = 0.5 * (Tw[faceI] + Trad_eff_local);
        scalar hrad_local = 4.0 * emissivity_ * sigma * pow3(Tmean);
        scalar heff_local = uValue_ + hrad_local;
        scalar Tref_local = (uValue_ * tempExt_ + hrad_local * Trad_eff_local + qSolar_[faceI] + qRadInt_[faceI]) / heff_local;
        scalar q_BC = heff_local * (Tw[faceI] - Tref_local);
        
        // The balance should be between OpenFOAM's calculated flux and the BC flux
        scalar balance_inner = qConvInt_[faceI] + qRadInt_[faceI] - q_BC;
        
        // For the boundary condition, we care about the balance at the inner surface
        qTotal_[faceI] = balance_inner;
    }
}


void Foam::buildingElementFvPatchScalarField::writeDebugInfo() const
{
    if (!debugTemp_ && !debugFlux_) return;
    
    const scalar currentTime = db().time().value();
    
    if ((currentTime - lastDebugTime_) < debugInterval_) return;
    
    lastDebugTime_ = currentTime;
    
    const scalarField Ti = this->patchInternalField();
    
    if (debugFace_ >= 0 && debugFace_ < Ti.size())
    {
        label i = debugFace_;
        
        Info<< nl << "===== BuildingElement Debug (" << patch().name() 
            << ", Face " << i << ") =====" << nl
            << "Time: " << currentTime << " s" << nl;
            
        if (debugTemp_)
        {
            // Get the cell temperature (not the patch internal field)
            const fvPatchField<scalar>& Tp = *this;
            const volScalarField& T = db().lookupObject<volScalarField>(internalField().name());
            const labelUList& faceCells = patch().faceCells();
            
            Info<< nl << "TEMPERATURES [C]:" << nl
                << "  T_cell:          " << T[faceCells[i]]-273.15 << nl
                << "  T_wall_inner:    " << Tp[i]-273.15 << nl
                << "  T_external:      " << tempExt_-273.15 << nl
                << "  T_sky:           " << tempSky_-273.15 << nl;
        }
            
        if (debugFlux_)
        {
            // For validation: also show the simple U-value based heat flux
            scalar q_validation = uValue_ * ((*this)[i] - tempExt_);
            
            Info<< nl << "HEAT FLUXES [W/m2]:" << nl
                << "  Interior side:" << nl
                << "    q_conv_int:    " << qConvInt_[i] 
                << " (>0: room->wall)" << nl
                << "    q_rad_int:     " << qRadInt_[i] 
                << " (net from OpenFOAM)" << nl
                << "  Through wall:" << nl
                << "    q_conduction:  " << qCond_[i] 
                << " (>0: inside->outside)" << nl
                << "    q_U-value:     " << q_validation
                << " (U × ΔT for validation)" << nl
                << "  Exterior side:" << nl
                << "    q_conv_ext:    " << qConvExt_[i] 
                << " (>0: air->wall)" << nl
                << "    q_rad_ext:     " << qRadExt_[i] 
                << " (>0: sky/env->wall)" << nl
                << "    q_solar:       " << qSolar_[i] 
                << " (>0: sun->wall)" << nl
                << "  Balance:" << nl
                << "    q_total:       " << qTotal_[i] 
                << " (should be ~0 at steady state)" << nl;
        }
        
        Info<< "==========================================" << nl << endl;
    }
    
    // Patch statistics
    if (debugFlux_)
    {
        Info<< "Patch statistics for " << patch().name() << ":" << nl
            << "  Area:            " << gSum(patch().magSf()) << " m2" << nl
            << "  Q_conduction:    " << gSum(qCond_ * patch().magSf()) << " W" << nl
            << "  Q_rad_internal:  " << gSum(qRadInt_ * patch().magSf()) 
            << " W (net)" << nl
            << "  Q_rad_external:  " << gSum(qRadExt_ * patch().magSf()) 
            << " W (net)" << nl
            << "  Q_solar:         " << gSum(qSolar_ * patch().magSf()) << " W" << nl
            << "  Q_total:         " << gSum(qTotal_ * patch().magSf()) << " W" << nl 
            << endl;
    }
}


void Foam::buildingElementFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalarField Ti = this->patchInternalField();
    
    if (useMultiLayer_)
    {
        // Update temperature distribution in layers
        updateLayerTemperatures();
        
        // Calculate heat flux based on temperature gradient at inner surface
        forAll(Ti, faceI)
        {
            const layer& firstLayer = layers_[0];
            const scalarField& Tfirst = Tlayers_[0];
            
            scalar T0 = Tfirst[faceI * firstLayer.nNodes];
            scalar T1 = Tfirst[faceI * firstLayer.nNodes + 1];
            scalar dx = firstLayer.thickness / (firstLayer.nNodes - 1);
            
            scalar gradT = (T1 - T0) / dx;
            qCond_[faceI] = -firstLayer.lambda * gradT;
        }
    }
    
    // Calculate all heat fluxes
    calculateHeatFluxes();
    
    // Write debug information if enabled
    writeDebugInfo();
    
    // Get thermal diffusivity at the boundary
    // For compressible flows, we need to calculate kappa from thermo properties
    scalarField kappa(patch().size());
    
    // Try different approaches depending on the solver
    if (db().foundObject<volScalarField>("alphaEff"))
    {
        const fvPatchScalarField& alphaEff = 
            patch().lookupPatchField<volScalarField, scalar>("alphaEff");
        kappa = alphaEff;
    }
    else if (db().foundObject<volScalarField>("kappaEff"))
    {
        const fvPatchScalarField& kappaEff = 
            patch().lookupPatchField<volScalarField, scalar>("kappaEff");
        kappa = kappaEff;
    }
    else
    {
        // For buoyantSimpleFoam and similar solvers, calculate from thermo
        // kappa = alpha * rho * Cp where alpha is thermal diffusivity
        // But for the BC, we can use a simplified approach
        
        // Get a reference value for thermal conductivity of air
        // k_air ≈ 0.026 W/(m·K) at room temperature
        scalar k_air = 0.026;
        
        // Get the cell spacing normal to the wall (not needed with our approach)
        
        // For the Robin BC formulation, we need kappa/delta
        // where kappa is thermal diffusivity (not conductivity)
        // Using kappa ≈ k_air as approximation
        kappa = k_air;
        
        Info<< "buildingElement BC: Using approximate thermal conductivity k = " 
            << k_air << " W/(m·K)" << endl;
    }
    
    // Solar gains (absorbed solar radiation)
    scalar qSolarGain = gValue_ * solarExt_;
    
    // Get internal radiation if available (ALREADY NET from OpenFOAM!)
    scalarField qRadiation(Ti.size(), 0.0);
    if (radField_ != "none" && db().foundObject<volScalarField>(radField_))
    {
        const volScalarField& qr = db().lookupObject<volScalarField>(radField_);
        qRadiation = qr.boundaryField()[patch().index()];
    }
    
    // External radiation and convection
    // The boundary condition formulation:
    // -k*dT/dn = U*(T_wall - T_ref) + q_sources
    // where q_sources includes solar gains and internal radiation
    const scalar sigma = physicoChemical::sigma.value();
    
    // Calculate boundary condition coefficients
    const scalarField& deltaCoeffs = patch().deltaCoeffs();
    
    // Mixed boundary condition formulation:
    // valueFraction*T + (1-valueFraction)*gradT = valueFraction*refValue + (1-valueFraction)*refGrad
    //
    // Our energy balance at the wall:
    // -k*dT/dn = h_eff*(T_wall - T_ref) - q_sources
    //
    // Where:
    // - h_eff includes U-value and radiation
    // - T_ref is the effective external temperature
    // - q_sources includes solar gains and internal radiation
    
    scalarField heff(Ti.size());
    scalarField Tref(Ti.size());
    
    forAll(Ti, faceI)
    {
        scalar Twall = this->operator[](faceI);
        scalar Trad_eff = viewSky_ * tempSky_ + (1 - viewSky_) * tempExt_;
        
        // Linearized radiation coefficient
        scalar Tmean = 0.5 * (Twall + Trad_eff);
        scalar hrad = 4.0 * emissivity_ * sigma * pow3(Tmean);
        
        // Total heat transfer coefficient (U-value + radiation)
        heff[faceI] = uValue_ + hrad;
        
        // Calculate reference temperature including all effects
        // Energy balance at wall: -k*dT/dn = h_eff*(T_wall - T_ref)
        // Where h_eff*(T_wall - T_ref) must equal all heat fluxes:
        // h_eff*(T_wall - T_ref) = U*(T_wall - T_ext) + h_rad*(T_wall - T_rad) - q_solar - q_rad_internal
        // 
        // Note: h_rad*(T_wall - T_rad) is positive when wall is warmer than radiation temperature
        // Note: q_solar is positive (heat gain), so we subtract it
        // Note: q_rad_internal from OpenFOAM is already net flux (positive into wall)
        //
        // Solving for T_ref:
        // h_eff*T_wall - h_eff*T_ref = U*T_wall - U*T_ext + h_rad*T_wall - h_rad*T_rad - q_solar - q_rad_internal
        // h_eff*T_ref = h_eff*T_wall - U*T_wall + U*T_ext - h_rad*T_wall + h_rad*T_rad + q_solar + q_rad_internal
        // h_eff*T_ref = T_wall*(h_eff - U - h_rad) + U*T_ext + h_rad*T_rad + q_solar + q_rad_internal
        // 
        // Since h_eff = U + h_rad:
        // h_eff*T_ref = T_wall*0 + U*T_ext + h_rad*T_rad + q_solar + q_rad_internal
        // 
        // Therefore:
        Tref[faceI] = (uValue_ * tempExt_ + hrad * Trad_eff + qSolarGain + qRadiation[faceI]) / heff[faceI];
    }
    
    // Set boundary condition coefficients
    valueFraction() = heff / (heff + kappa*deltaCoeffs);
    refValue() = Tref;
    refGrad() = 0.0;
    
    // Debug output for validation
    if (debugTemp_ && Pstream::master())
    {
        Info<< "BuildingElement BC Debug (updateCoeffs):" << nl
            << "  U-value: " << uValue_ << " W/(m²K)" << nl
            << "  kappa[0]: " << kappa[0] << " W/(m·K)" << nl
            << "  deltaCoeffs[0]: " << deltaCoeffs[0] << " 1/m" << nl
            << "  heff[0]: " << heff[0] << " W/(m²K)" << nl
            << "  valueFraction[0]: " << valueFraction()[0] << nl
            << "  refValue[0]: " << refValue()[0] << " K (" << refValue()[0] - 273.15 << " °C)" << nl
            << "  T_ext: " << tempExt_ << " K (" << tempExt_ - 273.15 << " °C)" << nl
            << endl;
    }
    
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::buildingElementFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    
    if (useMultiLayer_)
    {
        os.writeKeyword("layers") << nl;
        os << token::BEGIN_LIST << nl;
        forAll(layers_, i)
        {
            os << indent << token::BEGIN_BLOCK << nl;
            os.writeEntry("lambda", layers_[i].lambda);
            os.writeEntry("thickness", layers_[i].thickness);
            os.writeEntry("rho", layers_[i].rho);
            os.writeEntry("cp", layers_[i].cp);
            os.writeEntry("nNodes", layers_[i].nNodes);
            os << indent << token::END_BLOCK;
            if (i < layers_.size() - 1) os << token::COMMA;
            os << nl;
        }
        os << token::END_LIST << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeEntry("uValue", uValue_);
    }
    
    os.writeEntry("gValue", gValue_);
    os.writeEntry("Rsi", Rsi_);
    os.writeEntry("Rse", Rse_);
    os.writeEntry("tempExt", tempExt_);
    os.writeEntry("tempSky", tempSky_);
    os.writeEntry("solarExt", solarExt_);
    os.writeEntry("hExt", hExt_);
    os.writeEntry("emissivity", emissivity_);
    os.writeEntry("viewSky", viewSky_);
    os.writeEntry("radField", radField_);
    os.writeEntry("debugTemp", debugTemp_);
    os.writeEntry("debugFlux", debugFlux_);
    os.writeEntry("debugFace", debugFace_);
    os.writeEntry("debugInterval", debugInterval_);
    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        buildingElementFvPatchScalarField
    );
}

// ************************************************************************* //