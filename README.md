# HVAC-for-OpenFOAM
HVACTools for OpenFOAM - Make OpenFOAM for HVACTool Simulation great again.

Humidity Libary:
•	Switch to the following folder :
cd $FOAM_SRC/thermophysicalModels/basic

•	Add the following lines after the liquidThermo.C in under Make folder in files:
humidityRhoThermo/humidityRhoThermo.C
humidityRhoThermo/humidityRhoThermos.C 
humidityRhoThermo/derivedFvPatchFields/fixedHumidity/fixedHumidityFvPatchScalarField.C

•	Compile the thermo library and boundary condition 
wmake libso
