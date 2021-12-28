"comfortFoam" is OpenFOAM tool for thermal comfort according to DIN EN ISO 7730.

It is using a dictionary:
// 1 met = 58.15 W/m²	

clo 0.3; // Clothing [clo]
met	1.2; // Metabolic rate [met]
wme	0;	 // External work, normally around 0 [met]	
RH	50; // Relative humidity [%]

