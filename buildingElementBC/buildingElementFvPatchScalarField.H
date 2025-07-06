/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2412
    \\  /    A nd           | Web:      www.OpenFOAM.org
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef buildingElementFvPatchScalarField_H
#define buildingElementFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class buildingElementFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class buildingElementFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:
    // Layer structure for multi-layer walls
    struct layer
    {
        scalar lambda;      // Thermal conductivity [W/(mK)]
        scalar thickness;   // Layer thickness [m]
        scalar rho;         // Density [kg/m³]
        scalar cp;          // Specific heat capacity [J/(kgK)]
        label nNodes;       // Number of nodes for discretization
    };

private:
    // Private data

        //- Mode selection
        bool useMultiLayer_;        // Use multi-layer model instead of U-value
        
        //- Simple mode parameters
        scalar uValue_;             // Thermal transmittance [W/(m²K)]
        
        //- Multi-layer mode parameters
        List<layer> layers_;        // Layer properties
        PtrList<scalarField> Tlayers_;  // Temperature distribution in layers
        
        //- Common parameters
        scalar gValue_;             // Solar heat gain coefficient [-]
        scalar Rsi_;                // Internal surface resistance [m²K/W] - INFO ONLY (OpenFOAM calculates interior heat transfer)
        scalar Rse_;                // External surface resistance [m²K/W] - included in exterior calculation
        
        //- External boundary conditions
        scalar tempExt_;            // External air temperature [K]
        scalar tempSky_;            // Sky temperature [K]
        scalar solarExt_;           // External solar radiation [W/m²]
        scalar hExt_;               // External convection coefficient [W/(m²K)]
        scalar emissivity_;         // Surface emissivity [-]
        scalar viewSky_;            // View factor to sky [-]
        
        //- Internal radiation field name
        word radField_;
        
        //- Debug options
        bool debugTemp_;            // Debug temperature output
        bool debugFlux_;            // Debug heat flux output
        label debugFace_;           // Face index to monitor
        scalar debugInterval_;      // Output interval [s]
        mutable scalar lastDebugTime_;
        
        //- Heat flux fields for output
        mutable scalarField qConvInt_;
        mutable scalarField qRadInt_;
        mutable scalarField qCond_;
        mutable scalarField qConvExt_;
        mutable scalarField qRadExt_;
        mutable scalarField qSolar_;
        mutable scalarField qTotal_;


    // Private Member Functions

        //- Calculate effective U-value from layers
        scalar calculateUValue() const;
        
        //- Update multi-layer temperatures
        void updateLayerTemperatures();
        
        //- Calculate heat fluxes
        void calculateHeatFluxes() const;
        
        //- Write debug information
        void writeDebugInfo() const;


public:

    //- Runtime type information
    TypeName("buildingElement");


    // Constructors

        //- Construct from patch and internal field
        buildingElementFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        buildingElementFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given buildingElementFvPatchScalarField
        //  onto a new patch
        buildingElementFvPatchScalarField
        (
            const buildingElementFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        buildingElementFvPatchScalarField
        (
            const buildingElementFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new buildingElementFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        buildingElementFvPatchScalarField
        (
            const buildingElementFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new buildingElementFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //