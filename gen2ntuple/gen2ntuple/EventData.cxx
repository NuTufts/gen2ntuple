#include "EventData.h"
#include <algorithm>
#include <cstring>

namespace gen2ntuple {

// Static member definitions
constexpr int EventData::MAX_TRUE_PRIMARY_PARTS;
constexpr int EventData::MAX_TRUE_SIM_PARTS;
constexpr int EventData::MAX_KEYPOINTS;
constexpr int EventData::MAX_TRACKS;
constexpr int EventData::MAX_SHOWERS;

void EventData::clear() {
    // Event Identification
    fileid = -1;
    run = -1;
    subrun = -1;
    event = -1;
    
    // MC Truth Information
    xsecWeight = 1.0;
    
    // Neutrino truth
    trueNuE = -999.0;
    trueNuPDG = 0;
    trueNuCCNC = -1;
    trueNuMode = -1;
    trueNuIntrxnType = -1;
    
    // True lepton
    trueLepE = -999.0;
    trueLepPDG = 0;

    // True vertex
    trueVtxX = -999.0;
    trueVtxY = -999.0;
    trueVtxZ = -999.0;
    
    // True primary particles
    nTruePrimParts = 0;
    truePrimPartPDG.fill(0);

    truePrimPartE.fill(-999.0);
    truePrimPartPx.fill(-999.0);
    truePrimPartPy.fill(-999.0);
    truePrimPartPz.fill(-999.0);

    truePrimPartX.fill(-999.0);
    truePrimPartY.fill(-999.0);
    truePrimPartZ.fill(-999.0);

    truePrimPartContained.fill(-1);
    
    // True simulated particles
    nTrueSimParts = 0;
    trueSimPartPDG.fill(0);
    trueSimPartTID.fill(-1);
    trueSimPartMID.fill(-1);
    trueSimPartProcess.fill(-1);

    trueSimPartE.fill(-999.0);
    trueSimPartPx.fill(-999.0);
    trueSimPartPy.fill(-999.0);
    trueSimPartPz.fill(-999.0);

    trueSimPartX.fill(-999.0);
    trueSimPartY.fill(-999.0);
    trueSimPartZ.fill(-999.0);

    trueSimPartEDepX.fill(-999.0);
    trueSimPartEDepY.fill(-999.0);
    trueSimPartEDepZ.fill(-999.0);

    trueSimPartEndE.fill(-999.0);
    trueSimPartEndPx.fill(-999.0);
    trueSimPartEndPy.fill(-999.0);
    trueSimPartEndPz.fill(-999.0);

    trueSimPartEndX.fill(-999.0);
    trueSimPartEndY.fill(-999.0);
    trueSimPartEndZ.fill(-999.0);

    trueSimPartContained.fill(-1);

    trueSimPartPixelSumUplane.fill(0);
    trueSimPartPixelSumVplane.fill(0);
    trueSimPartPixelSumYplane.fill(0);

    
    // Vertex Reconstruction
    foundVertex = false;
    vtxIndex = -1;
    vtxScore = -999.0;
    vtxKPtype = -1;
    vtxKPscore = -999.0;
    
    // Vertex position
    vtxX = -999.0;
    vtxY = -999.0;
    vtxZ = -999.0;
    
    // Vertex quality
    vtxIsFiducial = false;
    vtxContainment = -999.0;
    vtxDistToTrue = -999.0;
    
    // Vertex metrics
    vtxMaxIntimePixelSum = -999.0;
    vtxFracHitsOnCosmic = -999.0;
    
    // Pixel fractions
    fracUnrecoIntimePixels.fill(-999.0);
    fracRecoOuttimePixels.fill(-999.0);
    
    // Keypoint Information
    nKeypoints = 0;
    kpClusterType.fill(-1);
    kpFilterType.fill(-1);
    kpMaxScore.fill(-999.0);
    kpMaxPosX.fill(-999.0);
    kpMaxPosY.fill(-999.0);
    kpMaxPosZ.fill(-999.0);
    
    // Track Information
    nTracks = 0;
    
    // Basic info
    trackIsSecondary.fill(false);
    trackNHits.fill(0);
    
    // Charge info
    trackCharge.fill(-999.0);
    trackChargeFrac.fill(-999.0);
    trackHitFrac.fill(-999.0);
    
    // Geometry
    trackStartPosX.fill(-999.0);
    trackStartPosY.fill(-999.0);
    trackStartPosZ.fill(-999.0);
    trackEndPosX.fill(-999.0);
    trackEndPosY.fill(-999.0);
    trackEndPosZ.fill(-999.0);
    trackStartDirX.fill(-999.0);
    trackStartDirY.fill(-999.0);
    trackStartDirZ.fill(-999.0);
    
    // Angular info
    trackCosTheta.fill(-999.0);
    trackCosThetaY.fill(-999.0);
    trackDistToVtx.fill(-999.0);
    
    // CNN scores
    trackElScore.fill(-999.0);
    trackPhScore.fill(-999.0);
    trackMuScore.fill(-999.0);
    trackPiScore.fill(-999.0);
    trackPrScore.fill(-999.0);
    trackPID.fill(0);
    
    // CNN quality
    trackComp.fill(-999.0);
    trackPurity.fill(-999.0);
    trackProcess.fill(-1);
    
    // Origin scores
    trackPrimaryScore.fill(-999.0);
    trackFromNeutralScore.fill(-999.0);
    trackFromChargedScore.fill(-999.0);
    
    // Energy
    trackRecoE.fill(-999.0);
    
    // Truth matching
    trackTruePDG.fill(0);
    trackTrueTID.fill(-1);
    trackTrueMID.fill(-1);
    trackTrueE.fill(-999.0);
    trackTrueComp.fill(-999.0);
    trackTruePurity.fill(-999.0);
    
    // Shower Information
    nShowers = 0;
    
    // Basic info
    showerIsSecondary.fill(false);
    showerNHits.fill(0);
    
    // Charge info
    showerCharge.fill(-999.0);
    showerChargeFrac.fill(-999.0);
    showerHitFrac.fill(-999.0);
    
    // Geometry
    showerStartPosX.fill(-999.0);
    showerStartPosY.fill(-999.0);
    showerStartPosZ.fill(-999.0);
    showerStartDirX.fill(-999.0);
    showerStartDirY.fill(-999.0);
    showerStartDirZ.fill(-999.0);
    
    // Angular info
    showerCosTheta.fill(-999.0);
    showerCosThetaY.fill(-999.0);
    showerDistToVtx.fill(-999.0);
    
    // CNN scores
    showerElScore.fill(-999.0);
    showerPhScore.fill(-999.0);
    showerMuScore.fill(-999.0);
    showerPiScore.fill(-999.0);
    showerPrScore.fill(-999.0);
    showerPID.fill(0);
    
    // CNN quality
    showerComp.fill(-999.0);
    showerPurity.fill(-999.0);
    showerProcess.fill(-1);
    
    // Origin scores
    showerPrimaryScore.fill(-999.0);
    showerFromNeutralScore.fill(-999.0);
    showerFromChargedScore.fill(-999.0);
    
    // Energy
    showerRecoE.fill(-999.0);
    showerRecoEU.fill(-999.0);
    showerRecoEV.fill(-999.0);
    showerRecoEY.fill(-999.0);
    
    // Truth matching
    showerTruePDG.fill(0);
    showerTrueTID.fill(-1);
    showerTrueMID.fill(-1);
    showerTrueE.fill(-999.0);
    showerTrueComp.fill(-999.0);
    showerTruePurity.fill(-999.0);
    
    // Event-level Features
    recoNuE = -999.0;
    
    // PCA analysis
    eventPCAxis0.fill(-999.0);
    eventPCAxis1.fill(-999.0);
    eventPCAxis2.fill(-999.0);
    eventPCEigenVals.fill(-999.0);
    eventPCProjMaxGap = -999.0;
    eventPCProjMaxDist = -999.0;
}

bool EventData::isMC() const {
    // Simple heuristic: if we have MC truth info, it's MC
    return (trueNuPDG != 0 || nTruePrimParts > 0 || nTrueSimParts > 0);
}

void POTData::clear() {
    totPOT = 0.0;
    totGoodPOT = 0.0;
}

} // namespace gen2ntuple