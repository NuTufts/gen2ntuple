#pragma once

#include <vector>
#include <array>

namespace gen2ntuple {

/**
 * @brief POD structure containing all event-level variables
 * 
 * This structure contains all the variables that will be saved to the ntuple,
 * organized by category. Each member corresponds to a branch in the output tree.
 */
struct EventData {
    
    // ===============================
    // Event Identification
    // ===============================
    int fileid;
    int run;
    int subrun; 
    int event;
    
    // ===============================
    // MC Truth Information (MC only)
    // ===============================

    // Xsec re-weight: filled by MCTruthProcessor::loadWeights
    float xsecWeight;
    
    // Neutrino truth: filled by MCTruthProcessor::extractNeutrinoTruth
    float trueNuE;
    int trueNuPDG;
    int trueNuCCNC;
    int trueNuMode;
    int trueNuIntrxnType;

    // True lepton: filled by MCTruthProcessor::extractNeutrinoTruth
    float trueLepE;
    int trueLepPDG;
    
    // True vertex: MCTruthProcessor::extractTrueVertex
    float trueVtxX;
    float trueVtxY;
    float trueVtxZ;
    
    // True primary particles (from GENIE): MCTruthProcessor::extractPrimaryParticles
    static constexpr int MAX_TRUE_PRIMARY_PARTS = 50;
    int nTruePrimParts;
    std::array<int, MAX_TRUE_PRIMARY_PARTS>   truePrimPartPDG;

    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartE;
    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartPx;
    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartPy;
    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartPz;

    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartX;
    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartY;
    std::array<float, MAX_TRUE_PRIMARY_PARTS> truePrimPartZ;

    std::array<int,   MAX_TRUE_PRIMARY_PARTS> truePrimPartContained;
    
    // True simulated particles (from Geant4): MCTruthProcessor::extractSimulatedParticles
    static constexpr int MAX_TRUE_SIM_PARTS = 200;
    int nTrueSimParts;
    std::array<int, MAX_TRUE_SIM_PARTS>   trueSimPartPDG;
    std::array<int, MAX_TRUE_SIM_PARTS>   trueSimPartTID;
    std::array<int, MAX_TRUE_SIM_PARTS>   trueSimPartMID;
    std::array<int, MAX_TRUE_SIM_PARTS>   trueSimPartProcess;

    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartE;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartPx;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartPy;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartPz;

    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartX;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartY;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartZ;

    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEDepX;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEDepY;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEDepZ;

    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndE;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndPx;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndPy;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndPz;

    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndX;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndY;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartEndZ;

    std::array<int, MAX_TRUE_SIM_PARTS>   trueSimPartContained;

    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartPixelSumUplane;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartPixelSumVplane;
    std::array<float, MAX_TRUE_SIM_PARTS> trueSimPartPixelSumYplane;
    
    // ===============================
    // Vertex Reconstruction
    // ===============================
    // VertexSelector::findBestVertex
    bool foundVertex;
    int vtxIndex;
    
    // Vertex position: filled by VertexSelector::findBestVertex
    float vtxX;
    float vtxY;
    float vtxZ;
    float vtxScore;
    
    // VertexSelector::checkFiducialVolume
    bool vtxIsFiducial;

    // VertexSelector::calculateContainment
    float vtxContainment;

    // VertexSelector::calculateTruthDistance
    float vtxDistToTrue;  // MC only
    
    // Vertex metrics: VertexSelector::calculateVertexQuality
    int vtxKPtype;
    float vtxKPscore;
    float vtxMaxIntimePixelSum;
    float vtxFracHitsOnCosmic;

    // Predicted Flash PE for the nu candidate
    // filled by VertexSelector::calculateFlashPrediction
    float observedPEtotal;
    std::array<float, 32> observedPE;

    float predictedPEtotal;
    std::array<float, 32> predictedPE;
    float sinkhorn_div;
    float fracerrPE;
    
    // Pixel fractions (one per plane)
    std::array<float, 3> fracUnrecoIntimePixels;
    std::array<float, 3> fracRecoOuttimePixels;
    
    // ===============================
    // Keypoint Information (optional)
    // ===============================
    static constexpr int MAX_KEYPOINTS = 100;
    int nKeypoints;
    std::array<int, MAX_KEYPOINTS> kpClusterType;
    std::array<int, MAX_KEYPOINTS> kpFilterType;
    std::array<float, MAX_KEYPOINTS> kpMaxScore;
    std::array<float, MAX_KEYPOINTS> kpMaxPosX;
    std::array<float, MAX_KEYPOINTS> kpMaxPosY;
    std::array<float, MAX_KEYPOINTS> kpMaxPosZ;
    
    // ===============================
    // Track Information
    // ===============================
    static constexpr int MAX_TRACKS = 50;
    int nTracks;
    
    // Basic info
    std::array<bool, MAX_TRACKS>  trackIsSecondary;
    std::array<int, MAX_TRACKS>   trackNHits;
    std::array<int, MAX_TRACKS>   trackIsContainedInFV;
    
    // Charge info: TrackProcessor::extractTrackInfo 
    std::array<float, MAX_TRACKS> trackCharge; //< TrackProcessor::calculateTrackCharge
    std::array<float, MAX_TRACKS> trackChargeFrac;
    std::array<float, MAX_TRACKS> trackHitFrac;
    
    // Geometry: filled by TrackProcessor::calculateTrackGeometry
    std::array<float, MAX_TRACKS> trackStartPosX;
    std::array<float, MAX_TRACKS> trackStartPosY;
    std::array<float, MAX_TRACKS> trackStartPosZ;
    std::array<float, MAX_TRACKS> trackEndPosX;
    std::array<float, MAX_TRACKS> trackEndPosY;
    std::array<float, MAX_TRACKS> trackEndPosZ;
    std::array<float, MAX_TRACKS> trackStartDirX;
    std::array<float, MAX_TRACKS> trackStartDirY;
    std::array<float, MAX_TRACKS> trackStartDirZ;
    
    // Angular info: filled by TrackProcessor::calculateTrackAngles
    std::array<float, MAX_TRACKS> trackCosTheta;
    std::array<float, MAX_TRACKS> trackCosThetaY;
    std::array<float, MAX_TRACKS> trackDistToVtx;

    // Classification status
    std::array<int, MAX_TRACKS> trackClassified;
    std::array<int, MAX_TRACKS> trackNGoodPlanes;
    
    // CNN scores
    std::array<float, MAX_TRACKS> trackElScore;
    std::array<float, MAX_TRACKS> trackPhScore;
    std::array<float, MAX_TRACKS> trackMuScore;
    std::array<float, MAX_TRACKS> trackPiScore;
    std::array<float, MAX_TRACKS> trackPrScore;
    std::array<int, MAX_TRACKS>   trackPID;
    
    // CNN quality
    std::array<float, MAX_TRACKS> trackComp;
    std::array<float, MAX_TRACKS> trackPurity;
    std::array<int, MAX_TRACKS>   trackProcess;
    
    // Origin scores
    std::array<float, MAX_TRACKS> trackPrimaryScore;
    std::array<float, MAX_TRACKS> trackFromNeutralScore;
    std::array<float, MAX_TRACKS> trackFromChargedScore;
    
    // Energy
    std::array<float, MAX_TRACKS> trackLength;
    std::array<float, MAX_TRACKS> trackRecoE;
    std::array<float, MAX_TRACKS> trackMuonE;
    std::array<float, MAX_TRACKS> trackProtonE;
    
    
    // Truth matching (MC only)
    std::array<int, MAX_TRACKS>   trackTruePDG;
    std::array<int, MAX_TRACKS>   trackTrueTID;
    std::array<int, MAX_TRACKS>   trackTrueMID;
    std::array<float, MAX_TRACKS> trackTrueE;
    std::array<float, MAX_TRACKS> trackTrueComp;
    std::array<float, MAX_TRACKS> trackTruePurity;
    
    // ===============================
    // Shower Information  
    // ===============================
    static constexpr int MAX_SHOWERS = 20;
    int nShowers;
    
    // Basic info
    std::array<bool, MAX_SHOWERS> showerIsSecondary;
    std::array<int, MAX_SHOWERS> showerNHits;
    std::array<int, MAX_SHOWERS> showerIsContainedInFV;
    
    // Charge info
    std::array<float, MAX_SHOWERS> showerCharge;
    std::array<float, MAX_SHOWERS> showerChargeFrac;
    std::array<float, MAX_SHOWERS> showerHitFrac;
    
    // Geometry
    std::array<float, MAX_SHOWERS> showerStartPosX;
    std::array<float, MAX_SHOWERS> showerStartPosY;
    std::array<float, MAX_SHOWERS> showerStartPosZ;
    std::array<float, MAX_SHOWERS> showerStartDirX;
    std::array<float, MAX_SHOWERS> showerStartDirY;
    std::array<float, MAX_SHOWERS> showerStartDirZ;
    
    // Angular info
    std::array<float, MAX_SHOWERS> showerCosTheta;
    std::array<float, MAX_SHOWERS> showerCosThetaY;
    std::array<float, MAX_SHOWERS> showerDistToVtx;
    
    // Classification status
    std::array<int, MAX_SHOWERS> showerClassified;
    std::array<int, MAX_SHOWERS> showerNGoodPlanes;
    
    // CNN scores
    std::array<float, MAX_SHOWERS> showerElScore;
    std::array<float, MAX_SHOWERS> showerPhScore;
    std::array<float, MAX_SHOWERS> showerMuScore;
    std::array<float, MAX_SHOWERS> showerPiScore;
    std::array<float, MAX_SHOWERS> showerPrScore;
    std::array<int, MAX_SHOWERS> showerPID;
    
    // CNN quality
    std::array<float, MAX_SHOWERS> showerComp;
    std::array<float, MAX_SHOWERS> showerPurity;
    std::array<int, MAX_SHOWERS> showerProcess;
    
    // Origin scores
    std::array<float, MAX_SHOWERS> showerPrimaryScore;
    std::array<float, MAX_SHOWERS> showerFromNeutralScore;
    std::array<float, MAX_SHOWERS> showerFromChargedScore;
    
    // Energy (multiple estimators)
    std::array<float, MAX_SHOWERS> showerRecoE;
    std::array<float, MAX_SHOWERS> showerRecoEU;
    std::array<float, MAX_SHOWERS> showerRecoEV;
    std::array<float, MAX_SHOWERS> showerRecoEY;
    std::array<float, MAX_SHOWERS> showerLength;
    std::array<float, MAX_SHOWERS> showerOpeningAngle;
    
    // Truth matching (MC only)
    std::array<int, MAX_SHOWERS> showerTruePDG;
    std::array<int, MAX_SHOWERS> showerTrueTID;
    std::array<int, MAX_SHOWERS> showerTrueMID;
    std::array<float, MAX_SHOWERS> showerTrueE;
    std::array<float, MAX_SHOWERS> showerTrueComp;
    std::array<float, MAX_SHOWERS> showerTruePurity;
    
    // ===============================
    // Event-level Features
    // ===============================
    float recoNuE;
    
    // PCA analysis
    std::array<float, 3> eventPCAxis0;
    std::array<float, 3> eventPCAxis1;
    std::array<float, 3> eventPCAxis2;
    std::array<float, 3> eventPCEigenVals;
    float eventPCProjMaxGap;
    float eventPCProjMaxDist;
    
    // ===============================
    // Initialization and utilities
    // ===============================
    
    /**
     * @brief Initialize all variables to default values
     */
    void clear();
    
    /**
     * @brief Check if this is MC data based on available truth info
     */
    bool isMC() const;
};

/**
 * @brief POD structure for POT information (separate tree for MC)
 */
struct POTData {
    double totPOT;
    double totGoodPOT;
    
    void clear();
};

} // namespace gen2ntuple