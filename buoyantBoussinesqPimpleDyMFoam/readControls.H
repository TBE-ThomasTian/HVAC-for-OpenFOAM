#include "readTimeControls.H"
#include "createControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", false);

checkMeshCourantNo = pimple.dict().lookupOrDefault("checkMeshCourantNo", false);
