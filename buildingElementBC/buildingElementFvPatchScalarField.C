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
    scalar Rtotal = 0.13;  // Internal surface resistance
    
    forAll(layers_, i)
    {
        Rtotal += layers_[i].thickness / layers_[i].lambda;
    }
    
    Rtotal += 0.04;  // External surface resistance
    
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
    const scalarField Ti = this->patchInternalField();
    const scalar sigma = physicoChemical::sigma.value();
    
    // Get internal radiation if available
    if (radField_ != "none" && db().foundObject<volScalarField>(radField_))
    {
        const volScalarField& qr = db().lookupObject<volScalarField>(radField_);
        // IMPORTANT: qr is already NET flux (incident - emitted) from OpenFOAM!
        qRadInt_ = qr.boundaryField()[patch().index()];
    }
    else
    {
        qRadInt_ = 0.0;
    }
    
    // Calculate heat fluxes for each face
    forAll(Ti, faceI)
    {
        // Wall inner temperature (this BC value)
        scalar Twall_int = this->operator[](faceI);
        
        // Outer surface temperature (simplified - in reality iterate)
        scalar Twall_ext = tempExt_;
        
        // Interior convection (simplified)
        scalar hint = 8.0; // W/(m²K) typical for interior
        qConvInt_[faceI] = hint * (Ti[faceI] - Twall_int);
        
        // Conduction through wall
        qCond_[faceI] = uValue_ * (Twall_int - tempExt_);
        
        // Exterior convection
        qConvExt_[faceI] = hExt_ * (Twall_ext - tempExt_);
        
        // Exterior radiation - WE calculate the NET flux
        scalar Trad_eff = viewSky_ * tempSky_ + (1 - viewSky_) * tempExt_;
        qRadExt_[faceI] = emissivity_ * sigma * (pow4(Twall_ext) - pow4(Trad_eff));
        
        // Solar gains
        qSolar_[faceI] = gValue_ * solarExt_;
        
        // Total balance
        // Note: qRadInt_ is already net (from OpenFOAM radiation model)
        qTotal_[faceI] = qConvInt_[faceI] + qRadInt_[faceI] 
                       - qCond_[faceI] - qRadExt_[faceI] + qSolar_[faceI];
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
            Info<< nl << "TEMPERATURES [C]:" << nl
                << "  T_fluid:         " << Ti[i]-273.15 << nl
                << "  T_wall_inner:    " << this->operator[](i)-273.15 << nl
                << "  T_external:      " << tempExt_-273.15 << nl
                << "  T_sky:           " << tempSky_-273.15 << nl;
        }
            
        if (debugFlux_)
        {
            Info<< nl << "HEAT FLUXES [W/m2] (positive = into wall):" << nl
                << "  q_conv_int:      " << qConvInt_[i] << nl
                << "  q_rad_int:       " << qRadInt_[i] 
                << " (net from OpenFOAM)" << nl
                << "  q_conduction:    " << qCond_[i] << nl
                << "  q_conv_ext:      " << qConvExt_[i] << nl
                << "  q_rad_ext:       " << qRadExt_[i] 
                << " (net to sky/environ)" << nl
                << "  q_solar:         " << qSolar_[i] << nl
                << "  q_total:         " << qTotal_[i] << nl;
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
    const fvPatchScalarField& kappa = 
        patch().lookupPatchField<volScalarField, scalar>("kappaEff");
    
    // Solar gains
    scalar qSolarGain = gValue_ * solarExt_;
    
    // Get internal radiation if available (ALREADY NET from OpenFOAM!)
    scalarField qRadiation(Ti.size(), 0.0);
    if (radField_ != "none" && db().foundObject<volScalarField>(radField_))
    {
        const volScalarField& qr = db().lookupObject<volScalarField>(radField_);
        qRadiation = qr.boundaryField()[patch().index()];
    }
    
    // External radiation (we calculate ourselves)
    const scalar sigma = physicoChemical::sigma.value();
    scalarField qRadExternal(Ti.size());
    forAll(Ti, faceI)
    {
        scalar Twall = this->operator[](faceI);
        scalar Trad_eff = viewSky_ * tempSky_ + (1 - viewSky_) * tempExt_;
        qRadExternal[faceI] = emissivity_ * sigma * (pow4(Twall) - pow4(Trad_eff));
    }
    
    // Robin boundary condition coefficients
    // k*dT/dn + alpha*T = beta
    scalarField alpha = uValue_ / kappa;
    scalarField beta = (uValue_*tempExt_ + qSolarGain + qRadiation - qRadExternal) / kappa;
    //                                                     ↑              ↑
    //                                              internal (net)   external (net)
    
    // Set boundary condition
    refValue() = beta/alpha;
    refGrad() = 0.0;
    valueFraction() = alpha/(alpha + patch().deltaCoeffs());
    
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