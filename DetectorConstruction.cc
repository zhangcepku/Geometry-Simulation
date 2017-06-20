#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include <G4SubtractionSolid.hh>
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"

#include "CLHEP/Units/SystemOfUnits.h"

DetectorConstruction::DetectorConstruction() {

// materials
//-----------
  DefineMaterials();
}


DetectorConstruction::~DetectorConstruction() {}


G4VPhysicalVolume* DetectorConstruction::Construct() {

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* trap_mat = nist->FindOrBuildMaterial("G4_Al");
  G4Material* pmt = nist->FindOrBuildMaterial("G4_C");
  

// Clean old geometry, if any
//----------------------------
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

// World
//=======

  G4Box*          solid  = new G4Box("Mother", 1000.0*CLHEP::cm, 1000.0*CLHEP::cm, 1000.0*CLHEP::cm);
  G4LogicalVolume*   logW   = new G4LogicalVolume(solid, pAir, "World");
  G4VPhysicalVolume* physW  = new G4PVPlacement(0, G4ThreeVector(), logW,
            "World", 0, false, 0);

  G4LogicalVolume*   logC   = new G4LogicalVolume(solid, pAir, "Calorimeter");
  new G4PVPlacement(0, G4ThreeVector(), logC, "Calorimeter",  logW, false, 0);

// PMT //
  G4double angle(4.5*CLHEP::deg), thick(2.0*CLHEP::cm), edge_1(35.3*CLHEP::cm), envthick(40.0*CLHEP::cm), anglep(45.*CLHEP::deg);

  G4double innerRadius = 0.*CLHEP::cm;
  G4double outerRadius = 2.1*CLHEP::cm;
  G4double hz = 19.3*CLHEP::cm;
  G4double trx = hz*cos(anglep);
  G4double startAngle = 0.*CLHEP::deg;
  G4double spanningAngle = 360.*CLHEP::deg;
  G4Tubs* solidcyl
    = new G4Tubs("PMT",
                 innerRadius,
                 outerRadius,
                 0.5*hz,
                 startAngle,
                 spanningAngle);

  G4LogicalVolume*   solidcyllog = new G4LogicalVolume(solidcyl, pmt, "PMT");
  G4double phiz_3 = 90.*CLHEP::deg + anglep - 0.5*angle;
  G4double phix_3 = phiz_3 + 90.0*CLHEP::deg;
  G4RotationMatrix* rot_3 = AddMatrix(90*CLHEP::deg, phix_3, 0, 0, 90*CLHEP::deg, phiz_3);
  G4double phiz_4 = 270.*CLHEP::deg + 45.*CLHEP::cm - 0.5*angle;
  G4double phix_4 = 360.*CLHEP::deg + 45.*CLHEP::cm - 0.5*angle;
  G4RotationMatrix* rot_4 = AddMatrix(90*CLHEP::deg, phix_4, 0, 0, 90*CLHEP::deg, phiz_4);

// env 1 //

  const unsigned int nc_1(5);
  std::string childNames_1[nc_1] = {"B8", "B9", "B10", "B11", "B12"};
  G4double    heights_1[nc_1]    = {66.8484225245044*CLHEP::cm, 
        37.8707804735234*CLHEP::cm, 42.3673111366067*CLHEP::cm, 
        55.8569031258564*CLHEP::cm, 
        64.8499644520229*CLHEP::cm};
  std::string env1Name("Envelope1");

  G4double toth_1(0);
  for (unsigned int k=0; k<nc_1; ++k) toth_1 += heights_1[k];
  G4double cfac = tan(0.5*angle);
  G4double bl1_1  = 0.5*edge_1;
  G4double bl2_1  = bl1_1 + toth_1*cfac;
  G4double h1_1  = 0.5*thick;

  G4Trap*  solid1 = new G4Trap(env1Name, 0.5*toth_1, 0, 0, 0.5*envthick, bl1_1, bl1_1, 0, 0.5*envthick, bl2_1, bl2_1, 0);
  G4LogicalVolume*   logE_1 = new G4LogicalVolume(solid1, pAir, env1Name);
  G4double zpos1 = -0.5*toth_1;
  G4double ypos1 = 0;
  for (unsigned int k=0; k<nc_1; ++k) {
    if (k==0) {
      ypos1 = 13*CLHEP::cm;
    } else if (k==1) {
      ypos1 = 3*CLHEP::cm;
    } else if (k==2) {
      ypos1 = 8*CLHEP::cm;
    } else if (k==3) {
      ypos1 = 13*CLHEP::cm;
    } else if (k==4) {
      ypos1 = 18*CLHEP::cm;
    }
    zpos1 += 0.5*heights_1[k];
    bl2_1   = bl1_1 + heights_1[k]*cfac;

    G4Trap*  outersolid1a = new G4Trap("outerchildNames_1[k]", 0.5*heights_1[k], 0, 0, h1_1, bl1_1, bl1_1, 0, h1_1, bl2_1, bl2_1, 0);
    G4Trap*  innersolid1a = new G4Trap("innerchildNames_1[k]", 0.5*heights_1[k] - 0.1*CLHEP::cm, 0, 0, h1_1 - 0.1*CLHEP::cm, bl1_1 - 0.1*CLHEP::cm, bl1_1 - 0.1*CLHEP::cm, 0, h1_1 - 0.1*CLHEP::cm, bl2_1 - 0.1*CLHEP::cm, bl2_1 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow1 = new G4SubtractionSolid("Hollow 1",outersolid1a,innersolid1a);
    G4LogicalVolume*   child_1 = new G4LogicalVolume(hollow1, trap_mat, "Hollow 1");

    new G4PVPlacement(0, G4ThreeVector(0, ypos1, zpos1), child_1, childNames_1[k],
          logE_1, false, 0);

    zpos1 += 0.5*heights_1[k];
    bl1_1   = bl2_1;
  }

  //Now the modules in mother
  G4double bl_1 = 0.5*tan(0.5*angle)*toth_1 + edge_1;
  G4double xpos_1 = 0.6*bl_1;
  G4double ypos_1 = -19*CLHEP::cm;
  G4double phiz_1 = 90.*CLHEP::deg-0.5*angle;
  G4double phix_1 = phiz_1 + 90.0*CLHEP::deg;
  G4RotationMatrix* rot_1 = AddMatrix(90*CLHEP::deg, phix_1, 0, 0, 90*CLHEP::deg, phiz_1);
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_1, ypos_1, -250.0*CLHEP::cm), logE_1, env1Name,
          logC, false, 0);

  // pmt //

  G4double ypos1p = -0.54*toth_1;
  G4double zpos1p = 0;
  G4double xpos1p = -6*CLHEP::cm;
  for (unsigned int k=0; k<nc_1; ++k) {
    if (k==0) {
      zpos1p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos1p = 3*CLHEP::cm;
    } else if (k==2) {
      zpos1p = 8*CLHEP::cm;
    } else if (k==3) {
      zpos1p = 13*CLHEP::cm;
    } else if (k==4) {
      zpos1p = 18*CLHEP::cm;
    }
    ypos1p += heights_1[k];

    new G4PVPlacement(rot_3, G4ThreeVector(xpos1p, ypos1p, zpos1p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_1x(0);
  const unsigned int nc_1x(5);
  std::string childNames_1x[nc_1x] = {"B8", "B9", "B10", "B11", "B12"};
  G4double    heights_1x[nc_1x]    =               {66.8484225245044*CLHEP::cm - 0.3*CLHEP::cm, 
        37.8707804735234*CLHEP::cm - 0.3*CLHEP::cm, 42.3673111366067*CLHEP::cm - 0.3*CLHEP::cm, 
        55.8569031258564*CLHEP::cm - 0.3*CLHEP::cm, 64.8499644520229*CLHEP::cm - 0.3*CLHEP::cm};
  std::string env1xName("Envelope1x");

  toth_1x = toth_1 - 0.3*CLHEP::cm;
  G4double bl1_1x  = 0.5*edge_1 - 0.15*CLHEP::cm;
  G4double bl2_1x  = bl1_1x + toth_1x*cfac;
  G4double h1_1x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid1x = new G4Trap(env1xName, 0.5*toth_1x, 0, 0, 0.5*envthick, bl1_1x, bl1_1x, 0, 0.5*envthick, bl2_1x, bl2_1x, 0);
  G4LogicalVolume*   logE_1x = new G4LogicalVolume(solid1x, pAir, env1xName);
  G4double zpos1x = -0.5*toth_1x;
  G4double ypos1x = 0;

  for (unsigned int k=0; k<nc_1x; ++k) {
    if (k==0) {
      ypos1x = 13*CLHEP::cm;
    } else if (k==1) {
      ypos1x = 3*CLHEP::cm;
    } else if (k==2) {
      ypos1x = 8*CLHEP::cm;
    } else if (k==3) {
      ypos1x = 13*CLHEP::cm;
    } else if (k==4) {
      ypos1x = 18*CLHEP::cm;
    }
    bl2_1x   = bl1_1x + heights_1x[k]*cfac;
    zpos1x += 0.5*heights_1x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid1xa = new G4Trap(childNames_1x[k], 0.5*heights_1x[k], 0, 0, h1_1x, bl1_1x, bl1_1x, 0, h1_1x, bl2_1x, bl2_1x, 0);
    G4LogicalVolume*   child_1x = new G4LogicalVolume(solid1xa, pSci, childNames_1x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos1x, zpos1x), child_1x, childNames_1x[k],
          logE_1x, false, 0);
    zpos1x += 0.5*heights_1x[k] + 0.15*CLHEP::cm;
    bl1_1x   = bl2_1x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_1x, env1xName,
          logE_1, false, 0);

// env 2 //

  const unsigned int nc_2(5);
  G4double edge_2(31*CLHEP::cm);
  std::string childNames_2[nc_2] = {"B7", "B8", "B9", "B10", "B11"};
  G4double    heights_2[nc_2]    = {55.5571344149842*CLHEP::cm, 
        66.8484225245044*CLHEP::cm, 37.8707804735234*CLHEP::cm, 
        42.3673111366067*CLHEP::cm, 55.8569031258564*CLHEP::cm};
  std::string env2Name("Envelope2");

  G4double toth_2(0);
  for (unsigned int k=0; k<nc_2; ++k) toth_2 += heights_2[k];
  G4double bl1_2  = 0.5*edge_2;
  G4double bl2_2  = bl1_2 + toth_2*cfac;
  G4double h1_2  = 0.5*thick;

  G4Trap* solid2 = new G4Trap(env2Name, 0.5*toth_2, 0, 0, 0.5*envthick, bl1_2, bl1_2, 0, 0.5*envthick, bl2_2, bl2_2, 0);
  G4LogicalVolume*   logE_2 = new G4LogicalVolume(solid2, pAir, env2Name);
  G4double zpos2 = -0.5*toth_2;
  G4double ypos2 = 0;
  for (unsigned int k=0; k<nc_2; ++k) {
    if (k==0) {
      ypos2 = 13*CLHEP::cm;
    } else if (k==1) {
      ypos2 = 8*CLHEP::cm;
    } else if (k==2) {
      ypos2 = 3*CLHEP::cm;
    } else if (k==3) {
      ypos2 = 18*CLHEP::cm;
    } else if (k==4) {
      ypos2 = 8*CLHEP::cm;
    }
    zpos2 += 0.5*heights_2[k];
    bl2_2   = bl1_2 + heights_2[k]*cfac;

    G4Trap*  outersolid2a = new G4Trap("outerchildNames_2[k]", 0.5*heights_2[k], 0, 0, h1_2, bl1_2, bl1_2, 0, h1_2, bl2_2, bl2_2, 0);
    G4Trap*  innersolid2a = new G4Trap("innerchildNames_2[k]", 0.5*heights_2[k] - 0.1*CLHEP::cm, 0, 0, h1_2 - 0.1*CLHEP::cm, bl1_2 - 0.1*CLHEP::cm, bl1_2 - 0.1*CLHEP::cm, 0, h1_2 - 0.1*CLHEP::cm, bl2_2 - 0.1*CLHEP::cm, bl2_2 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow2 = new G4SubtractionSolid("Hollow 2",outersolid2a,innersolid2a);
    G4LogicalVolume*   child_2 = new G4LogicalVolume(hollow2, trap_mat, "Hollow 2");

    new G4PVPlacement(0, G4ThreeVector(0, ypos2, zpos2), child_2, childNames_2[k],
          logE_2, false, 0);
    zpos2 += 0.5*heights_2[k];
    bl1_2   = bl2_2;
  }

  //Now the modules in mother
  G4double bl_2 = 0.5*tan(0.5*angle)*toth_2 + edge_2;
  G4double xpos_2 = 0.47*bl_2 + 2.1*xpos_1;
  G4double ypos_2 = ypos_1 + 0.4*(toth_1 - toth_2);
  G4double phiz_2 = 270.*CLHEP::deg-0.5*angle;
  G4double phix_2 = 360.*CLHEP::deg-0.5*angle;
  G4RotationMatrix* rot_2 = AddMatrix(90*CLHEP::deg, phix_2, 0, 0, 90*CLHEP::deg, phiz_2);
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_2, ypos_2, -250.0*CLHEP::cm), logE_2, env2Name,
          logC, false, 0);

  // pmt // 

  G4double ypos2p = 0.41*toth_2;
  G4double zpos2p = 0;
  G4double xpos2p = xpos1p + 2.15*bl_1 + trx;
  for (unsigned int k=0; k<nc_2; ++k) {
    if (k==0) {
      zpos2p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos2p = 8*CLHEP::cm;
    } else if (k==2) {
      zpos2p = 3*CLHEP::cm;
    } else if (k==3) {
      zpos2p = 18*CLHEP::cm;
    } else if (k==4) {
      zpos2p = 8*CLHEP::cm;
    }
    ypos2p += -heights_2[k];
    new G4PVPlacement(rot_4, G4ThreeVector(xpos2p, ypos2p, zpos2p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_2x(0);
  const unsigned int nc_2x(5);
  std::string childNames_2x[nc_2x] = {"B7", "B8", "B9", "B10", "B11"};
  G4double    heights_2x[nc_2x]    =               {55.5571344149842*CLHEP::cm - 0.3*CLHEP::cm, 
        66.8484225245044*CLHEP::cm - 0.3*CLHEP::cm, 37.8707804735234*CLHEP::cm - 0.3*CLHEP::cm, 
        42.3673111366067*CLHEP::cm - 0.3*CLHEP::cm, 55.8569031258564*CLHEP::cm - 0.3*CLHEP::cm};
  std::string env2xName("Envelope2x");

  toth_2x = toth_2 - 0.3*CLHEP::cm;
  G4double bl1_2x  = 0.5*edge_2 - 0.15*CLHEP::cm;
  G4double bl2_2x  = bl1_2x + toth_2x*cfac;
  G4double h1_2x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid2x = new G4Trap(env2xName, 0.5*toth_2x, 0, 0, 0.5*envthick, bl1_2x, bl1_2x, 0, 0.5*envthick, bl2_2x, bl2_2x, 0);
  G4LogicalVolume*   logE_2x = new G4LogicalVolume(solid2x, pAir, env2xName);
  G4double zpos2x = -0.5*toth_2x;
  G4double ypos2x = 0;

  for (unsigned int k=0; k<nc_2x; ++k) {
    if (k==0) {
      ypos2x = 13*CLHEP::cm;
    } else if (k==1) {
      ypos2x = 8*CLHEP::cm;
    } else if (k==2) {
      ypos2x = 3*CLHEP::cm;
    } else if (k==3) {
      ypos2x = 18*CLHEP::cm;
    } else if (k==4) {
      ypos2x = 8*CLHEP::cm;
    }
    bl2_2x   = bl1_2x + heights_2x[k]*cfac;
    zpos2x += 0.5*heights_2x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid2xa = new G4Trap(childNames_2x[k], 0.5*heights_2x[k], 0, 0, h1_2x, bl1_2x, bl1_2x, 0, h1_2x, bl2_2x, bl2_2x, 0);
    G4LogicalVolume*   child_2x = new G4LogicalVolume(solid2xa, pSci, childNames_2x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos2x, zpos2x), child_2x, childNames_2x[k],
          logE_2x, false, 0);
    zpos2x += 0.5*heights_2x[k] + 0.15*CLHEP::cm;
    bl1_2x   = bl2_2x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_2x, env2xName,
          logE_2, false, 0);

// env 3 //

  const unsigned int nc_3(4);
  G4double edge_3(42.8*CLHEP::cm);
  std::string childNames_3[nc_3] = {"C8", "C9", "C10", "C11"};
  G4double    heights_3[nc_3]    = {81.9367809717393*CLHEP::cm, 46.8638417996899*CLHEP::cm, 
        52.2596785953898*CLHEP::cm, 69.2465722114821*CLHEP::cm};
  std::string env3Name("Envelope3");

  G4double toth_3(0);
  for (unsigned int k=0; k<nc_3; ++k) toth_3 += heights_3[k];
  G4double bl1_3  = 0.5*edge_3;
  G4double bl2_3  = bl1_3 + toth_3*cfac;
  G4double h1_3  = 0.5*thick;

  G4Trap*  solid3 = new G4Trap(env3Name, 0.5*toth_3, 0, 0, 0.5*envthick, bl1_3, bl1_3, 0, 0.5*envthick, bl2_3, bl2_3, 0);
  G4LogicalVolume*   logE_3 = new G4LogicalVolume(solid3, pAir, env3Name);
  G4double zpos3 = -0.5*toth_3;
  G4double ypos3 = 0;
  for (unsigned int k=0; k<nc_3; ++k) {
    if (k==0) {
      ypos3 = 3*CLHEP::cm;
    } else if (k==1) {
      ypos3 = 18*CLHEP::cm;
    } else if (k==2) {
      ypos3 = 13*CLHEP::cm;
    } else if (k==3) {
      ypos3 = 18*CLHEP::cm;
    }
    zpos3 += 0.5*heights_3[k];
    bl2_3   = bl1_3 + heights_3[k]*cfac;

    G4Trap*  outersolid3a = new G4Trap("outerchildNames_3[k]", 0.5*heights_3[k], 0, 0, h1_3, bl1_3, bl1_3, 0, h1_3, bl2_3, bl2_3, 0);
    G4Trap*  innersolid3a = new G4Trap("innerchildNames_3[k]", 0.5*heights_3[k] - 0.1*CLHEP::cm, 0, 0, h1_3 - 0.1*CLHEP::cm, bl1_3 - 0.1*CLHEP::cm, bl1_3 - 0.1*CLHEP::cm, 0, h1_3 - 0.1*CLHEP::cm, bl2_3 - 0.1*CLHEP::cm, bl2_3 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow3 = new G4SubtractionSolid("Hollow 3",outersolid3a,innersolid3a);
    G4LogicalVolume*   child_3 = new G4LogicalVolume(hollow3, trap_mat, "Hollow 3");

    new G4PVPlacement(0, G4ThreeVector(0, ypos3, zpos3), child_3, childNames_3[k],
          logE_3, false, 0);
    zpos3 += 0.5*heights_3[k];
    bl1_3   = bl2_3;
  }

  //Now the modules in mother
  G4double bl_3 = 0.5*tan(0.5*angle)*toth_3 + edge_3;
  G4double xpos_3 = -0.53*bl_3;
  G4double ypos_3 = -10.086815 - 0.5*(toth_1 - toth_3);
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_3, ypos_3, -250.0*CLHEP::cm), logE_3, env3Name,
          logC, false, 0);

  // pmt //

  G4double ypos3p = 0.425*toth_3;
  G4double zpos3p = 0;
  G4double xpos3p = xpos1p + trx;
  for (unsigned int k=0; k<nc_3; ++k) {
    if (k==0) {
      zpos3p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos3p = 8*CLHEP::cm;
    } else if (k==2) {
      zpos3p = 3*CLHEP::cm;
    } else if (k==3) {
      zpos3p = 18*CLHEP::cm;
    } else if (k==4) {
      zpos3p = 8*CLHEP::cm;
    }
    ypos3p += -heights_3[k];
    new G4PVPlacement(rot_4, G4ThreeVector(xpos3p, ypos3p, zpos3p), solidcyllog, "pmt",
          logC, false, 0);
  }


  // Inside //

  G4double toth_3x(0);
  const unsigned int nc_3x(4);
  std::string childNames_3x[nc_3x] = {"C8", "C9", "C10", "C11"};
  G4double    heights_3x[nc_3x]    = {81.9367809717393*CLHEP::cm - 0.3*CLHEP::cm, 46.8638417996899*CLHEP::cm - 0.3*CLHEP::cm, 
        52.2596785953898*CLHEP::cm - 0.3*CLHEP::cm, 69.2465722114821*CLHEP::cm - 0.3*CLHEP::cm};
  std::string env3xName("Envelope3x");

  toth_3x = toth_3 - 0.3*CLHEP::cm;
  G4double bl1_3x  = 0.5*edge_3 - 0.15*CLHEP::cm;
  G4double bl2_3x  = bl1_3x + toth_3x*cfac;
  G4double h1_3x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid3x = new G4Trap(env3xName, 0.5*toth_3x, 0, 0, 0.5*envthick, bl1_3x, bl1_3x, 0, 0.5*envthick, bl2_3x, bl2_3x, 0);
  G4LogicalVolume*   logE_3x = new G4LogicalVolume(solid3x, pAir, env3xName);
  G4double zpos3x = -0.5*toth_3x;
  G4double ypos3x = 0;

  for (unsigned int k=0; k<nc_3x; ++k) {
    if (k==0) {
      ypos3x = 3*CLHEP::cm;
    } else if (k==1) {
      ypos3x = 18*CLHEP::cm;
    } else if (k==2) {
      ypos3x = 13*CLHEP::cm;
    } else if (k==3) {
      ypos3x = 18*CLHEP::cm;
    } 
    bl2_3x   = bl1_3x + heights_3x[k]*cfac;
    zpos3x += 0.5*heights_3x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid3xa = new G4Trap(childNames_3x[k], 0.5*heights_3x[k], 0, 0, h1_3x, bl1_3x, bl1_3x, 0, h1_3x, bl2_3x, bl2_3x, 0);
    G4LogicalVolume*   child_3x = new G4LogicalVolume(solid3xa, pSci, childNames_3x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos3x, zpos3x), child_3x, childNames_3x[k],
          logE_3x, false, 0);
    zpos3x += 0.5*heights_3x[k] + 0.15*CLHEP::cm;
    bl1_3x   = bl2_3x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_3x, env3xName,
          logE_3, false, 0);

// env 4 //

  const unsigned int nc_4(4);
  G4double edge_4 = edge_3;
  G4double toth_4 = toth_3;
  std::string childNames_4[nc_4] = {"C8", "C9", "C10", "C11"};
  G4double    heights_4[nc_4]    = {81.9367809717393*CLHEP::cm, 46.8638417996899*CLHEP::cm, 
        52.2596785953898*CLHEP::cm, 69.2465722114821*CLHEP::cm};
  std::string env4Name("Envelope4");

  G4double bl1_4  = 0.5*edge_4;
  G4double bl2_4  = bl1_4 + toth_4*cfac;
  G4double h1_4  = 0.5*thick;

  G4Trap*  solid4 = new G4Trap(env4Name, 0.5*toth_4, 0, 0, 0.5*envthick, bl1_4, bl1_4, 0, 0.5*envthick, bl2_4, bl2_4, 0);
  G4LogicalVolume*   logE_4 = new G4LogicalVolume(solid4, pAir, env4Name);
  G4double zpos4 = -0.5*toth_4;
  G4double ypos4 = 0;
  for (unsigned int k=0; k<nc_4; ++k) {
    if (k==0) {
      ypos4 = 3*CLHEP::cm;
    } else if (k==1) {
      ypos4 = 8*CLHEP::cm;
    } else if (k==2) {
      ypos4 = 13*CLHEP::cm;
    } else if (k==3) {
      ypos4 = 18*CLHEP::cm;
    }
    zpos4 += 0.5*heights_4[k];
    bl2_4   = bl1_4 + heights_4[k]*cfac;

    G4Trap*  outersolid4a = new G4Trap("outerchildNames_4[k]", 0.5*heights_4[k], 0, 0, h1_4, bl1_4, bl1_4, 0, h1_4, bl2_4, bl2_4, 0);
    G4Trap*  innersolid4a = new G4Trap("innerchildNames_4[k]", 0.5*heights_4[k] - 0.1*CLHEP::cm, 0, 0, h1_4 - 0.1*CLHEP::cm, bl1_4 - 0.1*CLHEP::cm, bl1_4 - 0.1*CLHEP::cm, 0, h1_4 - 0.1*CLHEP::cm, bl2_4 - 0.1*CLHEP::cm, bl2_4 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow4 = new G4SubtractionSolid("Hollow 4",outersolid4a,innersolid4a);
    G4LogicalVolume*   child_4 = new G4LogicalVolume(hollow4, trap_mat, "Hollow 4");

    new G4PVPlacement(0, G4ThreeVector(0, ypos4, zpos4), child_4, childNames_4[k],
          logE_4, false, 0);
    zpos4 += 0.5*heights_4[k];
    bl1_4   = bl2_4;
  }

  //Now the modules in mother
  G4double bl_4 = 0.5*tan(0.5*angle)*toth_4 + edge_4;
  G4double xpos_4 = -0.53*bl_4 + 2.1*xpos_3;
  G4double ypos_4 = ypos_3 + 2.0*CLHEP::cm;
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_4, ypos_4, -250.0*CLHEP::cm), logE_4, env4Name,
          logC, false, 0);

  // pmt //

  G4double ypos4p = -0.5*toth_4;
  G4double zpos4p = 0;
  G4double xpos4p = xpos3p - 2.5*bl_3;
  for (unsigned int k=0; k<nc_4; ++k) {
    if (k==0) {
      zpos4p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos4p = 3*CLHEP::cm;
    } else if (k==2) {
      zpos4p = 8*CLHEP::cm;
    } else if (k==3) {
      zpos4p = 13*CLHEP::cm;
    } else if (k==4) {
      zpos4p = 18*CLHEP::cm;
    }
    ypos4p += heights_4[k];

    new G4PVPlacement(rot_3, G4ThreeVector(xpos4p, ypos4p, zpos4p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_4x(0);
  const unsigned int nc_4x(4);
  std::string childNames_4x[nc_4x] = {"C8", "C9", "C10", "C11"};
  G4double    heights_4x[nc_4x]    = {81.9367809717393*CLHEP::cm - 0.3*CLHEP::cm, 46.8638417996899*CLHEP::cm - 0.3*CLHEP::cm, 
        52.2596785953898*CLHEP::cm - 0.3*CLHEP::cm, 69.2465722114821*CLHEP::cm - 0.3*CLHEP::cm};
  std::string env4xName("Envelope4x");

  toth_4x = toth_4 - 0.3*CLHEP::cm;
  G4double bl1_4x  = 0.5*edge_4 - 0.15*CLHEP::cm;
  G4double bl2_4x  = bl1_4x + toth_4x*cfac;
  G4double h1_4x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid4x = new G4Trap(env4xName, 0.5*toth_4x, 0, 0, 0.5*envthick, bl1_4x, bl1_4x, 0, 0.5*envthick, bl2_4x, bl2_4x, 0);
  G4LogicalVolume*   logE_4x = new G4LogicalVolume(solid4x, pAir, env4xName);
  G4double zpos4x = -0.5*toth_4x;
  G4double ypos4x = 0;

  for (unsigned int k=0; k<nc_4x; ++k) {
    if (k==0) {
      ypos4x = 3*CLHEP::cm;
    } else if (k==1) {
      ypos4x = 8*CLHEP::cm;
    } else if (k==2) {
      ypos4x = 13*CLHEP::cm;
    } else if (k==3) {
      ypos4x = 18*CLHEP::cm;
    } 
    bl2_4x   = bl1_4x + heights_4x[k]*cfac;
    zpos4x += 0.5*heights_4x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid4xa = new G4Trap(childNames_4x[k], 0.5*heights_4x[k], 0, 0, h1_4x, bl1_4x, bl1_4x, 0, h1_4x, bl2_4x, bl2_4x, 0);
    G4LogicalVolume*   child_4x = new G4LogicalVolume(solid4xa, pSci, childNames_4x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos4x, zpos4x), child_4x, childNames_4x[k],
          logE_4x, false, 0);
    zpos4x += 0.5*heights_4x[k] + 0.15*CLHEP::cm;
    bl1_4x   = bl2_4x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_4x, env4xName,
          logE_4, false, 0);

// env 5 //

  const unsigned int nc_5(6);
  G4double edge_5(20.3*CLHEP::cm);
  std::string childNames_5[nc_5] = {"C2", "C3", "C4", "C5", "C6", "C7"};
  G4double    heights_5[nc_5]    = {32.4749436778235*CLHEP::cm, 
        36.9714743409067*CLHEP::cm, 42.6670798474789*CLHEP::cm, 
        49.5617601975399*CLHEP::cm, 57.8553611983379*CLHEP::cm,
        68.9468035006099*CLHEP::cm};
  std::string env5Name("Envelope5");

  G4double toth_5(0);
  for (unsigned int k=0; k<nc_5; ++k) toth_5 += heights_5[k];
  G4double bl1_5  = 0.5*edge_5;
  G4double bl2_5  = bl1_5 + toth_5*cfac;
  G4double h1_5  = 0.5*thick;

  G4Trap*  solid5 = new G4Trap(env5Name, 0.5*toth_5, 0, 0, 0.5*envthick, bl1_5, bl1_5, 0, 0.5*envthick, bl2_5, bl2_5, 0);
  G4LogicalVolume*   logE_5 = new G4LogicalVolume(solid5, pAir, env5Name);
  G4double zpos5 = -0.5*toth_5;
  G4double ypos5 = 0;
  for (unsigned int k=0; k<nc_5; ++k) {
    if (k==0) {
      ypos5 = 3*CLHEP::cm;
    } else if (k==1) {
      ypos5 = 13*CLHEP::cm;
    } else if (k==2) {
      ypos5 = 8*CLHEP::cm;
    } else if (k==3) {
      ypos5 = 13*CLHEP::cm;
    } else if (k==4) {
      ypos5 = 3*CLHEP::cm;
    } else if (k==5) {
      ypos5 = 8*CLHEP::cm;
    }
    zpos5 += 0.5*heights_5[k];
    bl2_5   = bl1_5 + heights_5[k]*cfac;

    G4Trap*  outersolid5a = new G4Trap("outerchildNames_5[k]", 0.5*heights_5[k], 0, 0, h1_5, bl1_5, bl1_5, 0, h1_5, bl2_5, bl2_5, 0);
    G4Trap*  innersolid5a = new G4Trap("innerchildNames_5[k]", 0.5*heights_5[k] - 0.1*CLHEP::cm, 0, 0, h1_5 - 0.1*CLHEP::cm, bl1_5 - 0.1*CLHEP::cm, bl1_5 - 0.1*CLHEP::cm, 0, h1_5 - 0.1*CLHEP::cm, bl2_5 - 0.1*CLHEP::cm, bl2_5 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow5 = new G4SubtractionSolid("Hollow 5",outersolid5a,innersolid5a);
    G4LogicalVolume*   child_5 = new G4LogicalVolume(hollow5, trap_mat, "Hollow 5");

    new G4PVPlacement(0, G4ThreeVector(0, ypos5, zpos5), child_5, childNames_5[k],
          logE_5, false, 0);
    zpos5 += 0.5*heights_5[k];
    bl1_5   = bl2_5;
  }

  //Now the modules in mother
  G4double bl_5 = 0.5*tan(0.5*angle)*toth_5 + edge_5;
  G4double xpos_5 = 0.5*bl_5 + 1.08*xpos_2 + 0.5*bl_2;
  G4double ypos_5 = 5*CLHEP::cm;
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_5, ypos_5, -250.0*CLHEP::cm), logE_5, env5Name,
          logC, false, 0);

  // pmt //

  G4double ypos5p = -0.455*toth_5;
  G4double zpos5p = 0;
  G4double xpos5p = xpos2p - trx;
  for (unsigned int k=0; k<nc_5; ++k) {
    if (k==0) {
      zpos5p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos5p = 3*CLHEP::cm;
    } else if (k==2) {
      zpos5p = 8*CLHEP::cm;
    } else if (k==3) {
      zpos5p = 13*CLHEP::cm;
    } else if (k==4) {
      zpos5p = 18*CLHEP::cm;
    }
    ypos5p += heights_5[k];

    new G4PVPlacement(rot_3, G4ThreeVector(xpos5p, ypos5p, zpos5p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_5x(0);
  const unsigned int nc_5x(6);
  std::string childNames_5x[nc_5x] = {"C2", "C3", "C4", "C5", "C6", "C7"};
  G4double    heights_5x[nc_5x]    =               {32.4749436778235*CLHEP::cm - 0.3*CLHEP::cm, 
        36.9714743409067*CLHEP::cm - 0.3*CLHEP::cm, 42.6670798474789*CLHEP::cm - 0.3*CLHEP::cm, 
        49.5617601975399*CLHEP::cm - 0.3*CLHEP::cm, 57.8553611983379*CLHEP::cm - 0.3*CLHEP::cm,
        68.9468035006099*CLHEP::cm - 0.3*CLHEP::cm};
  std::string env5xName("Envelope5x");

  toth_5x = toth_5 - 0.3*CLHEP::cm;
  G4double bl1_5x  = 0.5*edge_5 - 0.15*CLHEP::cm;
  G4double bl2_5x  = bl1_5x + toth_5x*cfac;
  G4double h1_5x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid5x = new G4Trap(env5xName, 0.5*toth_5x, 0, 0, 0.5*envthick, bl1_5x, bl1_5x, 0, 0.5*envthick, bl2_5x, bl2_5x, 0);
  G4LogicalVolume*   logE_5x = new G4LogicalVolume(solid5x, pAir, env5xName);
  G4double zpos5x = -0.5*toth_5x;
  G4double ypos5x = 0;

  for (unsigned int k=0; k<nc_5x; ++k) {
    if (k==0) {
      ypos5x = 3*CLHEP::cm;
    } else if (k==1) {
      ypos5x = 13*CLHEP::cm;
    } else if (k==2) {
      ypos5x = 8*CLHEP::cm;
    } else if (k==3) {
      ypos5x = 13*CLHEP::cm;
    } else if (k==4) {
      ypos5x = 3*CLHEP::cm;
    } else if (k==5) {
      ypos5x = 8*CLHEP::cm;
    }
    bl2_5x   = bl1_5x + heights_5x[k]*cfac;
    zpos5x += 0.5*heights_5x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid5xa = new G4Trap(childNames_5x[k], 0.5*heights_5x[k], 0, 0, h1_5x, bl1_5x, bl1_5x, 0, h1_5x, bl2_5x, bl2_5x, 0);
    G4LogicalVolume*   child_5x = new G4LogicalVolume(solid5xa, pSci, childNames_5x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos5x, zpos5x), child_5x, childNames_5x[k],
          logE_5x, false, 0);
    zpos5x += 0.5*heights_5x[k] + 0.15*CLHEP::cm;
    bl1_5x   = bl2_5x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_5x, env5xName,
          logE_5, false, 0);

  // env 6 //

  const unsigned int nc_6(6);
  G4double edge_6(20.3*CLHEP::cm);
  std::string childNames_6[nc_6] = {"C2", "C3", "C4", "C5", "C6", "C7"};
  G4double    heights_6[nc_6]    = {32.4749436778235*CLHEP::cm, 
        36.9714743409067*CLHEP::cm, 42.6670798474789*CLHEP::cm, 
        49.5617601975399*CLHEP::cm, 57.8553611983379*CLHEP::cm,
        68.9468035006099*CLHEP::cm};
  std::string env6Name("Envelope6");

  G4double toth_6(0);
  for (unsigned int k=0; k<nc_6; ++k) toth_6 += heights_6[k];
  G4double bl1_6  = 0.5*edge_6;
  G4double bl2_6  = bl1_6 + toth_6*cfac;
  G4double h1_6  = 0.5*thick;

  G4Trap*  solid6 = new G4Trap(env6Name, 0.5*toth_6, 0, 0, 0.5*envthick, bl1_6, bl1_6, 0, 0.5*envthick, bl2_6, bl2_6, 0);
  G4LogicalVolume*   logE_6 = new G4LogicalVolume(solid6, pAir, env6Name);
  G4double zpos6 = -0.5*toth_6;
  G4double ypos6 = 0;
  for (unsigned int k=0; k<nc_6; ++k) {
    if (k==0) {
      ypos6 = 3*CLHEP::cm;
    } else if (k==1) {
      ypos6 = 13*CLHEP::cm;
    } else if (k==2) {
      ypos6 = 18*CLHEP::cm;
    } else if (k==3) {
      ypos6 = 8*CLHEP::cm;
    } else if (k==4) {
      ypos6 = 3*CLHEP::cm;
    } else if (k==5) {
      ypos6 = 18*CLHEP::cm;
    }
    zpos6 += 0.5*heights_6[k];
    bl2_6   = bl1_6 + heights_6[k]*cfac;

    G4Trap*  outersolid6a = new G4Trap("outerchildNames_6[k]", 0.5*heights_6[k], 0, 0, h1_6, bl1_6, bl1_6, 0, h1_6, bl2_6, bl2_6, 0);
    G4Trap*  innersolid6a = new G4Trap("innerchildNames_6[k]", 0.5*heights_6[k] - 0.1*CLHEP::cm, 0, 0, h1_6 - 0.1*CLHEP::cm, bl1_6 - 0.1*CLHEP::cm, bl1_6 - 0.1*CLHEP::cm, 0, h1_6 - 0.1*CLHEP::cm, bl2_6 - 0.1*CLHEP::cm, bl2_6 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow6 = new G4SubtractionSolid("Hollow 6",outersolid6a,innersolid6a);
    G4LogicalVolume*   child_6 = new G4LogicalVolume(hollow6, trap_mat, "Hollow 6");

    new G4PVPlacement(0, G4ThreeVector(0, ypos6, zpos6), child_6, childNames_6[k],
          logE_6, false, 0);
    zpos6 += 0.5*heights_6[k];
    bl1_6   = bl2_6;
  }

  //Now the modules in mother
  G4double xpos_6 = xpos_5 + 1.2*bl_5;
  G4double ypos_6 = ypos_5;
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_6, ypos_6, -250.0*CLHEP::cm), logE_6, env6Name,
          logC, false, 0);

  // pmt //

  G4double ypos6p = 0.49*toth_6;
  G4double zpos6p = 0;
  G4double xpos6p = xpos5p + 2.45*bl_5 + trx;
  for (unsigned int k=0; k<nc_6; ++k) {
    if (k==0) {
      zpos6p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos6p = 8*CLHEP::cm;
    } else if (k==2) {
      zpos6p = 3*CLHEP::cm;
    } else if (k==3) {
      zpos6p = 18*CLHEP::cm;
    } else if (k==4) {
      zpos6p = 8*CLHEP::cm;
    }
    ypos6p += -heights_6[k];
    new G4PVPlacement(rot_4, G4ThreeVector(xpos6p, ypos6p, zpos6p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_6x(0);
  const unsigned int nc_6x(6);
  std::string childNames_6x[nc_6x] = {"C2", "C3", "C4", "C5", "C6", "C7"};
  G4double    heights_6x[nc_6x]    =               {32.4749436778235*CLHEP::cm - 0.3*CLHEP::cm, 
        36.9714743409067*CLHEP::cm - 0.3*CLHEP::cm, 42.6670798474789*CLHEP::cm - 0.3*CLHEP::cm, 
        49.5617601975399*CLHEP::cm - 0.3*CLHEP::cm, 57.8553611983379*CLHEP::cm - 0.3*CLHEP::cm,
        68.9468035006099*CLHEP::cm - 0.3*CLHEP::cm};
  std::string env6xName("Envelope6x");

  toth_6x = toth_6 - 0.3*CLHEP::cm;
  G4double bl1_6x  = 0.5*edge_6 - 0.15*CLHEP::cm;
  G4double bl2_6x  = bl1_6x + toth_6x*cfac;
  G4double h1_6x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid6x = new G4Trap(env6xName, 0.5*toth_6x, 0, 0, 0.5*envthick, bl1_6x, bl1_6x, 0, 0.5*envthick, bl2_6x, bl2_6x, 0);
  G4LogicalVolume*   logE_6x = new G4LogicalVolume(solid6x, pAir, env6xName);
  G4double zpos6x = -0.5*toth_6x;
  G4double ypos6x = 0;

  for (unsigned int k=0; k<nc_6x; ++k) {
    if (k==0) {
      ypos6x = 3*CLHEP::cm;
    } else if (k==1) {
      ypos6x = 13*CLHEP::cm;
    } else if (k==2) {
      ypos6x = 18*CLHEP::cm;
    } else if (k==3) {
      ypos6x = 8*CLHEP::cm;
    } else if (k==4) {
      ypos6x = 3*CLHEP::cm;
    } else if (k==5) {
      ypos6x = 18*CLHEP::cm;
    }
    bl2_6x   = bl1_6x + heights_6x[k]*cfac;
    zpos6x += 0.5*heights_6x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid6xa = new G4Trap(childNames_6x[k], 0.5*heights_6x[k], 0, 0, h1_6x, bl1_6x, bl1_6x, 0, h1_6x, bl2_6x, bl2_6x, 0);
    G4LogicalVolume*   child_6x = new G4LogicalVolume(solid6xa, pSci, childNames_6x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos6x, zpos6x), child_6x, childNames_6x[k],
          logE_6x, false, 0);
    zpos6x += 0.5*heights_6x[k] + 0.15*CLHEP::cm;
    bl1_6x   = bl2_6x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_6x, env6xName,
          logE_6, false, 0);

// env 7 //

  const unsigned int nc_7(6);
  G4double edge_7(15.1*CLHEP::cm);
  std::string childNames_7[nc_7] = {"B1", "B2", "B3", "B4", "B5", "B6"};
  G4double    heights_7[nc_7]    = {23.3819594480329*CLHEP::cm, 
        26.4795694603792*CLHEP::cm, 30.2766397980939*CLHEP::cm, 
        34.6732475575531*CLHEP::cm, 40.4687759677493*CLHEP::cm,
        46.963764703314*CLHEP::cm};
  std::string env7Name("Envelope7");

  G4double toth_7(0);
  for (unsigned int k=0; k<nc_7; ++k) toth_7 += heights_7[k];
  G4double bl1_7  = 0.5*edge_7;
  G4double bl2_7  = bl1_7 + toth_7*cfac;
  G4double h1_7  = 0.5*thick;

  G4Trap*  solid7 = new G4Trap(env7Name, 0.5*toth_7, 0, 0, 0.5*envthick, bl1_7, bl1_7, 0, 0.5*envthick, bl2_7, bl2_7, 0);
  G4LogicalVolume*   logE_7 = new G4LogicalVolume(solid7, pAir, env7Name);
  G4double zpos7 = -0.5*toth_7;
  G4double ypos7 = 0;
  for (unsigned int k=0; k<nc_7; ++k) {
    if (k==0) {
      ypos7 = 3*CLHEP::cm;
    } else if (k==1) {
      ypos7 = 8*CLHEP::cm;
    } else if (k==2) {
      ypos7 = 3*CLHEP::cm;
    } else if (k==3) {
      ypos7 = 8*CLHEP::cm;
    } else if (k==4) {
      ypos7 = 3*CLHEP::cm;
    } else if (k==5) {
      ypos7 = 18*CLHEP::cm;
    }
    zpos7 += 0.5*heights_7[k];
    bl2_7  = bl1_7 + heights_7[k]*cfac;

    G4Trap*  outersolid7a = new G4Trap("outerchildNames_7[k]", 0.5*heights_7[k], 0, 0, h1_7, bl1_7, bl1_7, 0, h1_7, bl2_7, bl2_7, 0);
    G4Trap*  innersolid7a = new G4Trap("innerchildNames_7[k]", 0.5*heights_7[k] - 0.1*CLHEP::cm, 0, 0, h1_7 - 0.1*CLHEP::cm, bl1_7 - 0.1*CLHEP::cm, bl1_7 - 0.1*CLHEP::cm, 0, h1_7 - 0.1*CLHEP::cm, bl2_7 - 0.1*CLHEP::cm, bl2_7 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow7 = new G4SubtractionSolid("Hollow 7",outersolid7a,innersolid7a);
    G4LogicalVolume*   child_7 = new G4LogicalVolume(hollow7, trap_mat, "Hollow 7");

    new G4PVPlacement(0, G4ThreeVector(0, ypos7, zpos7), child_7, childNames_7[k],
          logE_7, false, 0);
    zpos7 += 0.5*heights_7[k];
    bl1_7   = bl2_7;
  }

  //Now the modules in mother
  G4double bl_7 = 0.5*tan(0.5*angle)*toth_7 + edge_7;
  G4double xpos_7 = -0.4*bl_7 + 4.3*xpos_3;
  G4double ypos_7 = ypos_3 + 0.54*50.0613747156602*CLHEP::cm;
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_7, ypos_7, -250.0*CLHEP::cm), logE_7, env7Name,
          logC, false, 0);

  // pmt //

  G4double ypos7p = 0.55*toth_7;
  G4double zpos7p = 0;
  G4double xpos7p = xpos4p + trx;
  for (unsigned int k=0; k<nc_7; ++k) {
    if (k==0) {
      zpos7p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos7p = 8*CLHEP::cm;
    } else if (k==2) {
      zpos7p = 3*CLHEP::cm;
    } else if (k==3) {
      zpos7p = 18*CLHEP::cm;
    } else if (k==4) {
      zpos7p = 8*CLHEP::cm;
    }
    ypos7p += -heights_7[k];
    new G4PVPlacement(rot_4, G4ThreeVector(xpos7p, ypos7p, zpos7p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_7x(0);
  const unsigned int nc_7x(6);
  std::string childNames_7x[nc_7x] = {"B1", "B2", "B3", "B4", "B5", "B6"};
  G4double    heights_7x[nc_7x]    =               {23.3819594480329*CLHEP::cm - 0.3*CLHEP::cm, 
        26.4795694603792*CLHEP::cm - 0.3*CLHEP::cm, 30.2766397980939*CLHEP::cm - 0.3*CLHEP::cm, 
        34.6732475575531*CLHEP::cm - 0.3*CLHEP::cm, 40.4687759677493*CLHEP::cm - 0.3*CLHEP::cm,
        46.963764703314*CLHEP::cm  - 0.3*CLHEP::cm};
  std::string env7xName("Envelope7x");

  toth_7x = toth_7 - 0.3*CLHEP::cm;
  G4double bl1_7x  = 0.5*edge_7 - 0.15*CLHEP::cm;
  G4double bl2_7x  = bl1_7x + toth_7x*cfac;
  G4double h1_7x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid7x = new G4Trap(env7xName, 0.5*toth_7x, 0, 0, 0.5*envthick, bl1_7x, bl1_7x, 0, 0.5*envthick, bl2_7x, bl2_7x, 0);
  G4LogicalVolume*   logE_7x = new G4LogicalVolume(solid7x, pAir, env7xName);
  G4double zpos7x = -0.5*toth_7x;
  G4double ypos7x = 0;

  for (unsigned int k=0; k<nc_7x; ++k) {
    if (k==0) {
      ypos7x = 3*CLHEP::cm;
    } else if (k==1) {
      ypos7x = 8*CLHEP::cm;
    } else if (k==2) {
      ypos7x = 3*CLHEP::cm;
    } else if (k==3) {
      ypos7x = 8*CLHEP::cm;
    } else if (k==4) {
      ypos7x = 3*CLHEP::cm;
    } else if (k==5) {
      ypos7x = 18*CLHEP::cm;
    }
    bl2_7x   = bl1_7x + heights_7x[k]*cfac;
    zpos7x += 0.5*heights_7x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid7xa = new G4Trap(childNames_7x[k], 0.5*heights_7x[k], 0, 0, h1_7x, bl1_7x, bl1_7x, 0, h1_7x, bl2_7x, bl2_7x, 0);
    G4LogicalVolume*   child_7x = new G4LogicalVolume(solid7xa, pSci, childNames_7x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos7x, zpos7x), child_7x, childNames_7x[k],
          logE_7x, false, 0);
    zpos7x += 0.5*heights_7x[k] + 0.15*CLHEP::cm;
    bl1_7x   = bl2_7x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_7x, env7xName,
          logE_7, false, 0);

// env 8 //

  const unsigned int nc_8(6);
  G4double edge_8(15.1*CLHEP::cm);
  std::string childNames_8[nc_8] = {"B1", "B2", "B3", "B4", "B5", "B6"};
  G4double    heights_8[nc_8]    = {23.3819594480329*CLHEP::cm, 
        26.4795694603792*CLHEP::cm, 30.2766397980939*CLHEP::cm, 
        34.6732475575531*CLHEP::cm, 40.4687759677493*CLHEP::cm,
        46.963764703314*CLHEP::cm};
  std::string env8Name("Envelope8");

  G4double toth_8(0);
  for (unsigned int k=0; k<nc_8; ++k) toth_8 += heights_8[k];
  G4double bl1_8  = 0.5*edge_8;
  G4double bl2_8  = bl1_8 + toth_8*cfac;
  G4double h1_8  = 0.5*thick;

  G4Trap*  solid8 = new G4Trap(env8Name, 0.5*toth_8, 0, 0, 0.5*envthick, bl1_8, bl1_8, 0, 0.5*envthick, bl2_8, bl2_8, 0);
  G4LogicalVolume*   logE_8 = new G4LogicalVolume(solid8, pAir, env8Name);
  G4double zpos8 = -0.5*toth_8;
  G4double ypos8 = 0;
  for (unsigned int k=0; k<nc_8; ++k) {
    if (k==0) {
      ypos8 = 3*CLHEP::cm;
    } else if (k==1) {
      ypos8 = 13*CLHEP::cm;
    } else if (k==2) {
      ypos8 = 8*CLHEP::cm;
    } else if (k==3) {
      ypos8 = 13*CLHEP::cm;
    } else if (k==4) {
      ypos8 = 18*CLHEP::cm;
    } else if (k==5) {
      ypos8 = 13*CLHEP::cm;
    }
    zpos8 += 0.5*heights_8[k];
    bl2_8   = bl1_8 + heights_8[k]*cfac;

    G4Trap*  outersolid8a = new G4Trap("outerchildNames_8[k]", 0.5*heights_8[k], 0, 0, h1_8, bl1_8, bl1_8, 0, h1_8, bl2_8, bl2_8, 0);
    G4Trap*  innersolid8a = new G4Trap("innerchildNames_8[k]", 0.5*heights_8[k] - 0.1*CLHEP::cm, 0, 0, h1_8 - 0.1*CLHEP::cm, bl1_8 - 0.1*CLHEP::cm, bl1_8 - 0.1*CLHEP::cm, 0, h1_8 - 0.1*CLHEP::cm, bl2_8 - 0.1*CLHEP::cm, bl2_8 - 0.1*CLHEP::cm, 0);
    G4SubtractionSolid *hollow8 = new G4SubtractionSolid("Hollow 8",outersolid8a,innersolid8a);
    G4LogicalVolume*   child_8 = new G4LogicalVolume(hollow8, trap_mat, "Hollow 8");

    new G4PVPlacement(0, G4ThreeVector(0, ypos8, zpos8), child_8, childNames_8[k],
          logE_8, false, 0);
    zpos8 += 0.5*heights_8[k];
    bl1_8   = bl2_8;
  }

  //Now the modules in mother
  G4double bl_8 = 0.5*tan(0.5*angle)*toth_8 + edge_8;
  G4double xpos_8 = -1.2*bl_8 + xpos_7;
  G4double ypos_8 = ypos_3 + 0.55*49.0613747156602*CLHEP::cm;
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_8, ypos_8, -250.0*CLHEP::cm), logE_8, env8Name,
          logC, false, 0);

  // pmt //

  G4double ypos8p = -0.38*toth_8;
  G4double zpos8p = 0;
  G4double xpos8p = xpos7p - 3.1*bl_7;
  for (unsigned int k=0; k<nc_8; ++k) {
    if (k==0) {
      zpos8p = 13*CLHEP::cm;
    } else if (k==1) {
      zpos8p = 3*CLHEP::cm;
    } else if (k==2) {
      zpos8p = 8*CLHEP::cm;
    } else if (k==3) {
      zpos8p = 13*CLHEP::cm;
    } else if (k==4) {
      zpos8p = 18*CLHEP::cm;
    }
    ypos8p += heights_8[k];

    new G4PVPlacement(rot_3, G4ThreeVector(xpos8p, ypos8p, zpos8p), solidcyllog, "pmt",
          logC, false, 0);
  }

  // Inside //

  G4double toth_8x(0);
  const unsigned int nc_8x(6);
  std::string childNames_8x[nc_8x] = {"B1", "B2", "B3", "B4", "B5", "B6"};
  G4double    heights_8x[nc_8x]    =               {23.3819594480329*CLHEP::cm - 0.3*CLHEP::cm, 
        26.4795694603792*CLHEP::cm - 0.3*CLHEP::cm, 30.2766397980939*CLHEP::cm - 0.3*CLHEP::cm, 
        34.6732475575531*CLHEP::cm - 0.3*CLHEP::cm, 40.4687759677493*CLHEP::cm - 0.3*CLHEP::cm,
        46.963764703314*CLHEP::cm  - 0.3*CLHEP::cm};
  std::string env8xName("Envelope8x");

  toth_8x = toth_8 - 0.3*CLHEP::cm;
  G4double bl1_8x  = 0.5*edge_8 - 0.15*CLHEP::cm;
  G4double bl2_8x  = bl1_8x + toth_8x*cfac;
  G4double h1_8x  = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid8x = new G4Trap(env8xName, 0.5*toth_8x, 0, 0, 0.5*envthick, bl1_8x, bl1_8x, 0, 0.5*envthick, bl2_8x, bl2_8x, 0);
  G4LogicalVolume*   logE_8x = new G4LogicalVolume(solid8x, pAir, env8xName);
  G4double zpos8x = -0.5*toth_8x;
  G4double ypos8x = 0;

  for (unsigned int k=0; k<nc_8x; ++k) {
    if (k==0) {
      ypos8x = 3*CLHEP::cm;
    } else if (k==1) {
      ypos8x = 13*CLHEP::cm;
    } else if (k==2) {
      ypos8x = 8*CLHEP::cm;
    } else if (k==3) {
      ypos8x = 13*CLHEP::cm;
    } else if (k==4) {
      ypos8x = 18*CLHEP::cm;
    } else if (k==5) {
      ypos8x = 13*CLHEP::cm;
    }
    bl2_8x  = bl1_8x + heights_8x[k]*cfac;
    zpos8x += 0.5*heights_8x[k] + 0.15*CLHEP::cm;

    G4Trap*  solid8xa = new G4Trap(childNames_8x[k], 0.5*heights_8x[k], 0, 0, h1_8x, bl1_8x, bl1_8x, 0, h1_8x, bl2_8x, bl2_8x, 0);
    G4LogicalVolume*   child_8x = new G4LogicalVolume(solid8xa, pSci, childNames_8x[k]);

    new G4PVPlacement(0, G4ThreeVector(0, ypos8x, zpos8x), child_8x, childNames_8x[k],
          logE_8x, false, 0);
    zpos8x += 0.5*heights_8x[k] + 0.15*CLHEP::cm;
    bl1_8x   = bl2_8x + 0.3*CLHEP::cm*cfac;
  }

  //Now the modules in mother
  new G4PVPlacement(0, G4ThreeVector(), logE_8x, env8xName,
          logE_8, false, 0);

// A9 //

  G4double edgeA9(25.8*CLHEP::cm);
  G4double bl1_A9 = 0.5*edgeA9;
  G4double bl2_A9 = 0.5*29.5*CLHEP::cm;
  G4double heightA9 = 50.0613747156602*CLHEP::cm;
  G4double h1_A9 = 0.5*thick; 

  G4Trap*  outersolidA9 = new G4Trap("outerA9", 0.5*heightA9, 0, 0, h1_A9, bl1_A9, bl1_A9, 0, h1_A9, bl2_A9, bl2_A9, 0);
  G4Trap*  innersolidA9 = new G4Trap("innerA9", 0.5*heightA9 - 0.1*CLHEP::cm, 0, 0, h1_A9 - 0.1*CLHEP::cm, bl1_A9 - 0.1*CLHEP::cm, bl1_A9 - 0.1*CLHEP::cm, 0, h1_A9 - 0.1*CLHEP::cm, bl2_A9 - 0.1*CLHEP::cm, bl2_A9 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollowA9 = new G4SubtractionSolid("Hollow A9",outersolidA9,innersolidA9);
  G4LogicalVolume*   solidlogA9 = new G4LogicalVolume(hollowA9, trap_mat, "A9");

  G4double xpos_A9 = 1.02*xpos_7;
  G4double ypos_A9 = -0.54*toth_7;
  
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_A9, ypos_A9, -237.0*CLHEP::cm), solidlogA9, "A9",
          logC, false, 0);


  G4double edgeA9x(25.8*CLHEP::cm);
  G4double bl1_A9x = 0.5*edgeA9x - 0.15*CLHEP::cm;
  G4double bl2_A9x = 0.5*29.5*CLHEP::cm - 0.15*CLHEP::cm;
  G4double heightA9x = 50.0613747156602*CLHEP::cm - 0.3*CLHEP::cm;
  G4double h1_A9x = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solidA9x = new G4Trap("A9x", 0.5*heightA9x, 0, 0, h1_A9x, bl1_A9x, bl1_A9x, 0, h1_A9x, bl2_A9x, bl2_A9x, 0);
  G4LogicalVolume*   solidlogA9x = new G4LogicalVolume(solidA9x, pSci, "A9x");
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_A9, ypos_A9, -237.0*CLHEP::cm), solidlogA9x, "A9x",
          logC, false, 0);

  G4double ypos_A9p = ypos7p - heightA9;
  G4double xpos_A9p = xpos7p;

  new G4PVPlacement(rot_4, G4ThreeVector(xpos_A9p, ypos_A9p, -237.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

// A8_1 //

  G4double edge_8_1(22.3*CLHEP::cm);
  G4double bl1_8_1 = 0.5*edge_8_1;
  G4double bl2_8_1 = 0.5*25.6*CLHEP::cm;
  G4double height8_1 = 40.9683904858696*CLHEP::cm;
  G4double h1_8_1 = 0.5*thick; 

  G4Trap*  outersolid8_1 = new G4Trap("outerA8_1", 0.5*height8_1, 0, 0, h1_8_1, bl1_8_1, bl1_8_1, 0, h1_8_1, bl2_8_1, bl2_8_1, 0);
  G4Trap*  innersolid8_1 = new G4Trap("innerA8_1", 0.5*height8_1 - 0.1*CLHEP::cm, 0, 0, h1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, 0, h1_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollow8_1 = new G4SubtractionSolid("Hollow A8_1",outersolid8_1,innersolid8_1);
  G4LogicalVolume*   solidlog8_1 = new G4LogicalVolume(hollow8_1, trap_mat, "A8_1");

  G4double xpos_8_1 = 1.125*xpos_2;
  G4double ypos_8_1 = 0.519*toth_2;
  
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_8_1, ypos_8_1, -247.0*CLHEP::cm), solidlog8_1, "A8_1",
          logC, false, 0);

  G4double edge_8_1x(22.3*CLHEP::cm);
  G4double bl1_8_1x = 0.5*edge_8_1x - 0.15*CLHEP::cm;
  G4double bl2_8_1x = 0.5*25.6*CLHEP::cm - 0.15*CLHEP::cm;
  G4double height8_1x = 40.9683904858696*CLHEP::cm - 0.3*CLHEP::cm;
  G4double h1_8_1x = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid8_1x = new G4Trap("A8_1x", 0.5*height8_1x, 0, 0, h1_8_1x, bl1_8_1x, bl1_8_1x, 0, h1_8_1x, bl2_8_1x, bl2_8_1x, 0);
  G4LogicalVolume*   solidlog8_1x = new G4LogicalVolume(solid8_1x, pSci, "A8_1x");
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_8_1, ypos_8_1, -247.0*CLHEP::cm), solidlog8_1x, "A8_1x",
          logC, false, 0);

  G4double xpos_8_1p = xpos2p;
  G4double ypos_8_1p = ypos_8_1 - 0.65*height8_1;

  new G4PVPlacement(rot_4, G4ThreeVector(xpos_8_1p, ypos_8_1p, -247.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

// A8_2 //

  G4Trap*  outersolid8_2 = new G4Trap("outerA8_2", 0.5*height8_1, 0, 0, h1_8_1, bl1_8_1, bl1_8_1, 0, h1_8_1, bl2_8_1, bl2_8_1, 0);
  G4Trap*  innersolid8_2 = new G4Trap("innerA8_2", 0.5*height8_1 - 0.1*CLHEP::cm, 0, 0, h1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, 0, h1_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollow8_2 = new G4SubtractionSolid("Hollow A8_2",outersolid8_2,innersolid8_2);
  G4LogicalVolume*   solidlog8_2 = new G4LogicalVolume(hollow8_2, trap_mat, "A8_2");

  G4double bl_8_2 = 0.5*tan(0.5*angle)*height8_1 + edge_8_1;
  G4double xpos_8_2 = xpos_8_1 - 1.05*bl_8_2;
  G4double ypos_8_2 = ypos_8_1;
  
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_8_2, ypos_8_2, -242.0*CLHEP::cm), solidlog8_2, "A8_2",
          logC, false, 0);

  G4Trap*  solid8_2x = new G4Trap("A8_2x", 0.5*height8_1x, 0, 0, h1_8_1x, bl1_8_1x, bl1_8_1x, 0, h1_8_1x, bl2_8_1x, bl2_8_1x, 0);
  G4LogicalVolume*   solidlog8_2x = new G4LogicalVolume(solid8_2x, pSci, "A8_2x");
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_8_2, ypos_8_2, -242.0*CLHEP::cm), solidlog8_2x, "A8_2x",
          logC, false, 0);

  G4double xpos_8_2p = xpos_8_1p - 2.1*bl_8_2 - trx;
  G4double ypos_8_2p = ypos_8_1p + height8_1 + trx;

  new G4PVPlacement(rot_3, G4ThreeVector(xpos_8_2p, ypos_8_2p, -242.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

// A8_3 //

  G4Trap*  outersolid8_3 = new G4Trap("outerA8_3", 0.5*height8_1, 0, 0, h1_8_1, bl1_8_1, bl1_8_1, 0, h1_8_1, bl2_8_1, bl2_8_1, 0);
  G4Trap*  innersolid8_3 = new G4Trap("innerA8_3", 0.5*height8_1 - 0.1*CLHEP::cm, 0, 0, h1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, 0, h1_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollow8_3 = new G4SubtractionSolid("Hollow A8_3",outersolid8_3,innersolid8_3);
  G4LogicalVolume*   solidlog8_3 = new G4LogicalVolume(hollow8_3, trap_mat, "A8_3");

  G4double xpos_8_3 = xpos_8_2 - 1.05*bl_8_2;
  G4double ypos_8_3 = ypos_8_1 + 1*CLHEP::cm;
  
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_8_3, ypos_8_3, -247.0*CLHEP::cm), solidlog8_3, "A8_3",
          logC, false, 0);

  G4Trap*  solid8_3x = new G4Trap("A8_3x", 0.5*height8_1x, 0, 0, h1_8_1x, bl1_8_1x, bl1_8_1x, 0, h1_8_1x, bl2_8_1x, bl2_8_1x, 0);
  G4LogicalVolume*   solidlog8_3x = new G4LogicalVolume(solid8_3x, pSci, "A8_3x");
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_8_3, ypos_8_3, -247.0*CLHEP::cm), solidlog8_3x, "A8_3x",
          logC, false, 0);

  G4double ypos_8_3p = ypos_8_1p;
  G4double xpos_8_3p = xpos_8_2p + trx;

  new G4PVPlacement(rot_4, G4ThreeVector(xpos_8_3p, ypos_8_3p, -247.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

// A8_4 //

  G4Trap*  outersolid8_4 = new G4Trap("outerA8_4", 0.5*height8_1, 0, 0, h1_8_1, bl1_8_1, bl1_8_1, 0, h1_8_1, bl2_8_1, bl2_8_1, 0);
  G4Trap*  innersolid8_4 = new G4Trap("innerA8_4", 0.5*height8_1 - 0.1*CLHEP::cm, 0, 0, h1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, bl1_8_1 - 0.1*CLHEP::cm, 0, h1_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, bl2_8_1 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollow8_4 = new G4SubtractionSolid("Hollow A8_4",outersolid8_4,innersolid8_4);
  G4LogicalVolume*   solidlog8_4 = new G4LogicalVolume(hollow8_4, trap_mat, "A8_4");

  G4double xpos_8_4 = xpos_8_3 - 1.05*bl_8_2;
  G4double ypos_8_4 = 1.01*ypos_8_1;
  
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_8_4, ypos_8_4, -237.0*CLHEP::cm), solidlog8_4, "A8_4",
          logC, false, 0);

  G4Trap*  solid8_4x = new G4Trap("A8_4x", 0.5*height8_1x, 0, 0, h1_8_1x, bl1_8_1x, bl1_8_1x, 0, h1_8_1x, bl2_8_1x, bl2_8_1x, 0);
  G4LogicalVolume*   solidlog8_4x = new G4LogicalVolume(solid8_4x, pSci, "A8_4x");
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_8_4, ypos_8_4, -237.0*CLHEP::cm), solidlog8_4x, "A8_4x",
          logC, false, 0);

  G4double ypos_8_4p = ypos_8_2p;
  G4double xpos_8_4p = xpos_8_3p - 2*bl_8_2 - trx;

  new G4PVPlacement(rot_3, G4ThreeVector(xpos_8_4p, ypos_8_4p, -237.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

// A7_1 //

  G4double edge_7_1(19.8*CLHEP::cm);
  G4double bl1_7_1 = 0.5*edge_7_1;
  G4double bl2_7_1 = 0.5*22.3*CLHEP::cm;
  G4double height7_1 = 34.1736330394327*CLHEP::cm;
  G4double h1_7_1 = 0.5*thick; 

  G4Trap*  outersolid7_1 = new G4Trap("outerA7_1", 0.5*height7_1, 0, 0, h1_7_1, bl1_7_1, bl1_7_1, 0, h1_7_1, bl2_7_1, bl2_7_1, 0);
  G4Trap*  innersolid7_1 = new G4Trap("innerA7_1", 0.5*height7_1 - 0.1*CLHEP::cm, 0, 0, h1_7_1 - 0.1*CLHEP::cm, bl1_7_1 - 0.1*CLHEP::cm, bl1_7_1 - 0.1*CLHEP::cm, 0, h1_7_1 - 0.1*CLHEP::cm, bl2_7_1 - 0.1*CLHEP::cm, bl2_7_1 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollow7_1 = new G4SubtractionSolid("Hollow A7_1",outersolid7_1,innersolid7_1);
  G4LogicalVolume*   solidlog7_1 = new G4LogicalVolume(hollow7_1, trap_mat, "A7_1");

  G4double xpos_7_1 = xpos_8_4 - 0.95*bl_8_2;
  G4double ypos_7_1 = 0.99*ypos_8_1;
  
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_7_1, ypos_7_1, -232.0*CLHEP::cm), solidlog7_1, "A7_1",
          logC, false, 0);

  G4double edge_7_1x(19.8*CLHEP::cm);
  G4double bl1_7_1x = 0.5*edge_7_1x - 0.15*CLHEP::cm;
  G4double bl2_7_1x = 0.5*22.3*CLHEP::cm - 0.15*CLHEP::cm;
  G4double height7_1x = 34.1736330394327*CLHEP::cm - 0.3*CLHEP::cm;
  G4double h1_7_1x = 0.5*thick - 0.15*CLHEP::cm;

  G4Trap*  solid7_1x = new G4Trap("A7_1x", 0.5*height7_1x, 0, 0, h1_7_1x, bl1_7_1x, bl1_7_1x, 0, h1_7_1x, bl2_7_1x, bl2_7_1x, 0);
  G4LogicalVolume*   solidlog7_1x = new G4LogicalVolume(solid7_1x, pSci, "A7_1x");
  new G4PVPlacement(rot_2, G4ThreeVector(xpos_7_1, ypos_7_1, -232.0*CLHEP::cm), solidlog7_1x, "A7_1x",
          logC, false, 0);

  G4double ypos_7_1p = ypos_7_1 - 0.75*height7_1;
  G4double xpos_7_1p = xpos_8_4p + 0.9*trx;

  new G4PVPlacement(rot_4, G4ThreeVector(xpos_7_1p, ypos_7_1p, -232.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

// A7_2 //

  G4Trap*  outersolid7_2 = new G4Trap("outerA7_2", 0.5*height7_1, 0, 0, h1_7_1, bl1_7_1, bl1_7_1, 0, h1_7_1, bl2_7_1, bl2_7_1, 0);
  G4Trap*  innersolid7_2 = new G4Trap("innerA7_2", 0.5*height7_1 - 0.1*CLHEP::cm, 0, 0, h1_7_1 - 0.1*CLHEP::cm, bl1_7_1 - 0.1*CLHEP::cm, bl1_7_1 - 0.1*CLHEP::cm, 0, h1_7_1 - 0.1*CLHEP::cm, bl2_7_1 - 0.1*CLHEP::cm, bl2_7_1 - 0.1*CLHEP::cm, 0);
  G4SubtractionSolid *hollow7_2 = new G4SubtractionSolid("Hollow A7_2",outersolid7_2,innersolid7_2);
  G4LogicalVolume*   solidlog7_2 = new G4LogicalVolume(hollow7_2, trap_mat, "A7_2");

  G4double bl_7_2 = 0.5*tan(0.5*angle)*height7_1 + edge_7_1;
  G4double xpos_7_2 = xpos_7_1 - bl_7_2;
  G4double ypos_7_2 = 0.99*ypos_8_1;
  
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_7_2, ypos_7_2, -242.0*CLHEP::cm), solidlog7_2, "A7_2",
          logC, false, 0);

  G4Trap*  solid7_2x = new G4Trap("A7_2x", 0.5*height7_1x, 0, 0, h1_7_1x, bl1_7_1x, bl1_7_1x, 0, h1_7_1x, bl2_7_1x, bl2_7_1x, 0);
  G4LogicalVolume*   solidlog7_2x = new G4LogicalVolume(solid7_2x, pSci, "A7_2x");
  new G4PVPlacement(rot_1, G4ThreeVector(xpos_7_2, ypos_7_2, -242.0*CLHEP::cm), solidlog7_2x, "A7_2x",
          logC, false, 0);

  G4double ypos_7_2p = ypos_8_4p - 0.7*(height8_1 - height7_1);
  G4double xpos_7_2p = xpos_7_1p - 2*bl_7_2 - trx;

  new G4PVPlacement(rot_3, G4ThreeVector(xpos_7_2p, ypos_7_2p, -242.0*CLHEP::cm), solidcyllog, "pmt",
          logC, false, 0);

  return physW;  
}


void DetectorConstruction::DefineMaterials() { 

  //
  // define Elements
  //-----------------

  G4Element* H  = new G4Element("Hydrogen","H", 1.,  1.01*CLHEP::g/CLHEP::mole);
  G4Element* C  = new G4Element("Carbon"  ,"C", 6., 12.01*CLHEP::g/CLHEP::mole);
  G4Element* N  = new G4Element("Nitrogen","N", 7., 14.01*CLHEP::g/CLHEP::mole);
  G4Element* O  = new G4Element("Oxygen"  ,"O", 8., 16.00*CLHEP::g/CLHEP::mole);

  // define scintillator (C_9H_10)_n
  //---------------------------------
  pSci = new G4Material("Scintillator", 1.032*CLHEP::g/CLHEP::cm3, 2);
  pSci->AddElement(C, 9);
  pSci->AddElement(H, 10);

  const G4int nSci = 1;
  G4double eSci[nSci] = { 3.10*CLHEP::eV };
  G4double rSci[nSci] = { 1.58    };
 
  G4MaterialPropertiesTable* proSci = new G4MaterialPropertiesTable();
  proSci->AddProperty("RINDEX", eSci, rSci, nSci);
  pSci->SetMaterialPropertiesTable(proSci);


  // define Air:
  //------------
  pAir = new G4Material("Air", 1.290*CLHEP::mg/CLHEP::cm3, 2);
  pAir->AddElement(N, 0.7);
  pAir->AddElement(O, 0.3);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


G4RotationMatrix* DetectorConstruction::AddMatrix(G4double th1, 
						  G4double phi1, 
						  G4double th2, 
						  G4double phi2, 
						  G4double th3, 
						  G4double phi3) {

  G4double sinth1 = std::sin(th1); 
  G4double costh1 = std::cos(th1);
  G4double sinth2 = std::sin(th2);
  G4double costh2 = std::cos(th2);
  G4double sinth3 = std::sin(th3);
  G4double costh3 = std::cos(th3);

  G4double sinph1 = std::sin(phi1); 
  G4double cosph1 = std::cos(phi1);
  G4double sinph2 = std::sin(phi2);
  G4double cosph2 = std::cos(phi2);
  G4double sinph3 = std::sin(phi3);
  G4double cosph3 = std::cos(phi3);
				    
  //xprime axis coordinates
  CLHEP::Hep3Vector xprime(sinth1*cosph1,sinth1*sinph1,costh1);
  //yprime axis coordinates
  CLHEP::Hep3Vector yprime(sinth2*cosph2,sinth2*sinph2,costh2);
  //zprime axis coordinates
  CLHEP::Hep3Vector zprime(sinth3*cosph3,sinth3*sinph3,costh3);

  G4RotationMatrix *rotMat = new G4RotationMatrix();
  rotMat->rotateAxes(xprime, yprime, zprime);
  if (*rotMat == G4RotationMatrix()) {
    delete rotMat;
    rotMat = 0;
  } else {
    rotMat->invert();
  }

  return rotMat;
}
