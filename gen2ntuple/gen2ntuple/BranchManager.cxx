#include "BranchManager.h"
#include <iostream>
#include <sstream>

namespace gen2ntuple {

BranchManager::BranchManager(TFile* output_file, bool is_mc, bool include_keypoints)
    : output_file_(output_file), is_mc_(is_mc), include_keypoints_(include_keypoints), 
      event_tree_(nullptr), pot_tree_(nullptr) {
    
    // Create the event tree (owned by TFile)
    event_tree_ = new TTree("EventTree", "Gen2 Analysis Ntuple");
    
    // Create POT tree if MC (owned by TFile)
    if (is_mc_) {
        pot_tree_ = new TTree("potTree", "POT Information");
    }
}

BranchManager::~BranchManager() = default;

void BranchManager::setupEventBranches(EventData* event_data) {
    if (!event_tree_ || !event_data) {
        std::cerr << "Error: Event tree or data pointer is null" << std::endl;
        return;
    }
    
    setupEventIDBranches(event_data);
    
    if (is_mc_) {
        setupMCTruthBranches(event_data);
    }
    
    setupVertexBranches(event_data);
    
    if (include_keypoints_) {
        setupKeypointBranches(event_data);
    }
    
    setupTrackBranches(event_data);
    setupShowerBranches(event_data);
    setupEventFeatureBranches(event_data);
}

void BranchManager::setupPOTBranches(POTData* pot_data) {
    if (!is_mc_ || !pot_tree_ || !pot_data) {
        if (!is_mc_) {
            std::cerr << "Warning: POT tree requested but not MC data" << std::endl;
        }
        return;
    }
    
    pot_tree_->Branch("totPOT", &pot_data->totPOT, "totPOT/D");
    pot_tree_->Branch("totGoodPOT", &pot_data->totGoodPOT, "totGoodPOT/D");
}

void BranchManager::setupEventIDBranches(EventData* data) {
    event_tree_->Branch("fileid", &data->fileid, "fileid/I");
    event_tree_->Branch("run", &data->run, "run/I");
    event_tree_->Branch("subrun", &data->subrun, "subrun/I");
    event_tree_->Branch("event", &data->event, "event/I");
}

void BranchManager::setupMCTruthBranches(EventData* data) {
    // Cross-section weight
    event_tree_->Branch("xsecWeight", &data->xsecWeight, "xsecWeight/F");
    
    // Neutrino truth
    event_tree_->Branch("trueNuE", &data->trueNuE, "trueNuE/F");
    event_tree_->Branch("trueNuPDG", &data->trueNuPDG, "trueNuPDG/I");
    event_tree_->Branch("trueNuCCNC", &data->trueNuCCNC, "trueNuCCNC/I");
    event_tree_->Branch("trueNuMode", &data->trueNuMode, "trueNuMode/I");
    event_tree_->Branch("trueNuIntrxnType", &data->trueNuIntrxnType, "trueNuIntrxnType/I");
    
    // True lepton
    event_tree_->Branch("trueLepE",   &data->trueLepE,   "trueLepE/F");
    event_tree_->Branch("trueLepPDG", &data->trueLepPDG, "trueLepPDG/I");

    // True vertex
    event_tree_->Branch("trueVtxX", &data->trueVtxX, "trueVtxX/F");
    event_tree_->Branch("trueVtxY", &data->trueVtxY, "trueVtxY/F");
    event_tree_->Branch("trueVtxZ", &data->trueVtxZ, "trueVtxZ/F");
    
    // True primary particles
    event_tree_->Branch("nTruePrimParts", &data->nTruePrimParts, "nTruePrimParts/I");
    createArrayBranch(event_tree_, "truePrimPartPDG", data->truePrimPartPDG.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "I");
    createArrayBranch(event_tree_, "truePrimPartE", data->truePrimPartE.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    createArrayBranch(event_tree_, "truePrimPartPx", data->truePrimPartPx.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    createArrayBranch(event_tree_, "truePrimPartPy", data->truePrimPartPy.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    createArrayBranch(event_tree_, "truePrimPartPz", data->truePrimPartPz.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    
    createArrayBranch(event_tree_, "truePrimPartX", data->truePrimPartX.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    createArrayBranch(event_tree_, "truePrimPartY", data->truePrimPartY.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    createArrayBranch(event_tree_, "truePrimPartZ", data->truePrimPartZ.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "F");
    createArrayBranch(event_tree_, "truePrimPartContained", data->truePrimPartContained.data(), 
                     "nTruePrimParts", EventData::MAX_TRUE_PRIMARY_PARTS, "I");
    
    // True simulated particles
    event_tree_->Branch("nTrueSimParts", &data->nTrueSimParts, "nTrueSimParts/I");
    createArrayBranch(event_tree_, "trueSimPartPDG", data->trueSimPartPDG.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "I");
    createArrayBranch(event_tree_, "trueSimPartTID", data->trueSimPartTID.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "I");
    createArrayBranch(event_tree_, "trueSimPartMID", data->trueSimPartMID.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "I");
    createArrayBranch(event_tree_, "trueSimPartProcess", data->trueSimPartProcess.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "I");
    
    createArrayBranch(event_tree_, "trueSimPartE", data->trueSimPartE.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartPx", data->trueSimPartPx.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartPy", data->trueSimPartPy.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartPz", data->trueSimPartPz.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");

    createArrayBranch(event_tree_, "trueSimPartX", data->trueSimPartX.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartY", data->trueSimPartY.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartZ", data->trueSimPartZ.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");


    createArrayBranch(event_tree_, "trueSimPartEDepX", data->trueSimPartEDepX.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEDepY", data->trueSimPartEDepY.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEDepZ", data->trueSimPartEDepZ.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");

    createArrayBranch(event_tree_, "trueSimPartEndE", data->trueSimPartEndE.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEndPx", data->trueSimPartEndPx.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEndPy", data->trueSimPartEndPy.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEndPz", data->trueSimPartEndPz.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");

    createArrayBranch(event_tree_, "trueSimPartEndX", data->trueSimPartEndX.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEndY", data->trueSimPartEndY.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartEndZ", data->trueSimPartEndZ.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");

    createArrayBranch(event_tree_, "trueSimPartContained", data->trueSimPartContained.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "I");

    createArrayBranch(event_tree_, "trueSimPartPixelSumUplane", data->trueSimPartPixelSumUplane.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartPixelSumVplane", data->trueSimPartPixelSumVplane.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");
    createArrayBranch(event_tree_, "trueSimPartPixelSumYplane", data->trueSimPartPixelSumYplane.data(), 
                     "nTrueSimParts", EventData::MAX_TRUE_SIM_PARTS, "F");

}

void BranchManager::setupVertexBranches(EventData* data) {
    // Vertex reconstruction
    event_tree_->Branch("foundVertex", &data->foundVertex, "foundVertex/O");
    event_tree_->Branch("vtxIndex", &data->vtxIndex, "vtxIndex/I");
    
    // Vertex position
    event_tree_->Branch("vtxX", &data->vtxX, "vtxX/F");
    event_tree_->Branch("vtxY", &data->vtxY, "vtxY/F");
    event_tree_->Branch("vtxZ", &data->vtxZ, "vtxZ/F");
    event_tree_->Branch("vtxScore", &data->vtxScore, "vtxScore/F");
    
    // Vertex quality
    event_tree_->Branch("vtxIsFiducial", &data->vtxIsFiducial, "vtxIsFiducial/O");
    event_tree_->Branch("vtxContainment", &data->vtxContainment, "vtxContainment/F");
    
    if (is_mc_) {
        event_tree_->Branch("vtxDistToTrue", &data->vtxDistToTrue, "vtxDistToTrue/F");
    }

    // Vertex metrics
    event_tree_->Branch("vtxKPtype", &data->vtxKPtype, "vtxKPtype/I");
    event_tree_->Branch("vtxKPscore", &data->vtxKPscore, "vtxKPscore/F");
    event_tree_->Branch("vtxMaxIntimePixelSum", &data->vtxMaxIntimePixelSum, "vtxMaxIntimePixelSum/F");
    event_tree_->Branch("vtxFracHitsOnCosmic", &data->vtxFracHitsOnCosmic, "vtxFracHitsOnCosmic/F");

    // Flash Predictions
    event_tree_->Branch("observedPEtotal", &data->observedPEtotal, "observedPEtotal/F");
    event_tree_->Branch("observedPE", data->observedPE.data(), "observedPE[32]/F");
    event_tree_->Branch("predictedPEtotal", &data->predictedPEtotal, "predictedPEtotal/F");
    event_tree_->Branch("predictedPE", data->predictedPE.data(), "predictedPE[32]/F" );
    event_tree_->Branch("sinkhorn_div", &data->sinkhorn_div, "sinkhorn_div/F");
    event_tree_->Branch("fracerrPE", &data->fracerrPE, "fracerrPE/F");

    
    // Pixel fractions (fixed size arrays)
    event_tree_->Branch("fracUnrecoIntimePixels", data->fracUnrecoIntimePixels.data(), "fracUnrecoIntimePixels[3]/F");
    event_tree_->Branch("fracRecoOuttimePixels", data->fracRecoOuttimePixels.data(), "fracRecoOuttimePixels[3]/F");
}

void BranchManager::setupKeypointBranches(EventData* data) {
    event_tree_->Branch("nKeypoints", &data->nKeypoints, "nKeypoints/I");
    createArrayBranch(event_tree_, "kpClusterType", data->kpClusterType.data(), 
                     "nKeypoints", EventData::MAX_KEYPOINTS, "I");
    createArrayBranch(event_tree_, "kpFilterType", data->kpFilterType.data(), 
                     "nKeypoints", EventData::MAX_KEYPOINTS, "I");
    createArrayBranch(event_tree_, "kpMaxScore", data->kpMaxScore.data(), 
                     "nKeypoints", EventData::MAX_KEYPOINTS, "F");
    createArrayBranch(event_tree_, "kpMaxPosX", data->kpMaxPosX.data(), 
                     "nKeypoints", EventData::MAX_KEYPOINTS, "F");
    createArrayBranch(event_tree_, "kpMaxPosY", data->kpMaxPosY.data(), 
                     "nKeypoints", EventData::MAX_KEYPOINTS, "F");
    createArrayBranch(event_tree_, "kpMaxPosZ", data->kpMaxPosZ.data(), 
                     "nKeypoints", EventData::MAX_KEYPOINTS, "F");
}

void BranchManager::setupTrackBranches(EventData* data) {
    event_tree_->Branch("nTracks", &data->nTracks, "nTracks/I");
    
    // Basic info
    createArrayBranch(event_tree_, "trackIsSecondary", data->trackIsSecondary.data(), 
                     "nTracks", EventData::MAX_TRACKS, "O");
    createArrayBranch(event_tree_, "trackNHits", data->trackNHits.data(), 
                     "nTracks", EventData::MAX_TRACKS, "I");
    
    // Charge info
    createArrayBranch(event_tree_, "trackCharge", data->trackCharge.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackChargeFrac", data->trackChargeFrac.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackHitFrac", data->trackHitFrac.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    
    // Geometry
    createArrayBranch(event_tree_, "trackStartPosX", data->trackStartPosX.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackStartPosY", data->trackStartPosY.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackStartPosZ", data->trackStartPosZ.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackEndPosX", data->trackEndPosX.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackEndPosY", data->trackEndPosY.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackEndPosZ", data->trackEndPosZ.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackStartDirX", data->trackStartDirX.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackStartDirY", data->trackStartDirY.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackStartDirZ", data->trackStartDirZ.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    
    // Angular info
    createArrayBranch(event_tree_, "trackCosTheta", data->trackCosTheta.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackCosThetaY", data->trackCosThetaY.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackDistToVtx", data->trackDistToVtx.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    
    // CNN scores
    createArrayBranch(event_tree_, "trackElScore", data->trackElScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackPhScore", data->trackPhScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackMuScore", data->trackMuScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackPiScore", data->trackPiScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackPrScore", data->trackPrScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackPID", data->trackPID.data(), 
                     "nTracks", EventData::MAX_TRACKS, "I");
    
    // CNN quality
    createArrayBranch(event_tree_, "trackComp", data->trackComp.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackPurity", data->trackPurity.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackProcess", data->trackProcess.data(), 
                     "nTracks", EventData::MAX_TRACKS, "I");
    
    // Origin scores
    createArrayBranch(event_tree_, "trackPrimaryScore", data->trackPrimaryScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackFromNeutralScore", data->trackFromNeutralScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    createArrayBranch(event_tree_, "trackFromChargedScore", data->trackFromChargedScore.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    
    // Energy
    createArrayBranch(event_tree_, "trackRecoE", data->trackRecoE.data(), 
                     "nTracks", EventData::MAX_TRACKS, "F");
    
    // Truth matching (MC only)
    if (is_mc_) {
        createArrayBranch(event_tree_, "trackTruePDG", data->trackTruePDG.data(), 
                         "nTracks", EventData::MAX_TRACKS, "I");
        createArrayBranch(event_tree_, "trackTrueTID", data->trackTrueTID.data(), 
                         "nTracks", EventData::MAX_TRACKS, "I");
        createArrayBranch(event_tree_, "trackTrueMID", data->trackTrueMID.data(), 
                         "nTracks", EventData::MAX_TRACKS, "I");
        createArrayBranch(event_tree_, "trackTrueE", data->trackTrueE.data(), 
                         "nTracks", EventData::MAX_TRACKS, "F");
        createArrayBranch(event_tree_, "trackTrueComp", data->trackTrueComp.data(), 
                         "nTracks", EventData::MAX_TRACKS, "F");
        createArrayBranch(event_tree_, "trackTruePurity", data->trackTruePurity.data(), 
                         "nTracks", EventData::MAX_TRACKS, "F");
    }
}

void BranchManager::setupShowerBranches(EventData* data) {
    event_tree_->Branch("nShowers", &data->nShowers, "nShowers/I");
    
    // Basic info
    createArrayBranch(event_tree_, "showerIsSecondary", data->showerIsSecondary.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "O");
    createArrayBranch(event_tree_, "showerNHits", data->showerNHits.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "I");
    
    // Charge info
    createArrayBranch(event_tree_, "showerCharge", data->showerCharge.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerChargeFrac", data->showerChargeFrac.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerHitFrac", data->showerHitFrac.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    
    // Geometry
    createArrayBranch(event_tree_, "showerStartPosX", data->showerStartPosX.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerStartPosY", data->showerStartPosY.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerStartPosZ", data->showerStartPosZ.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerStartDirX", data->showerStartDirX.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerStartDirY", data->showerStartDirY.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerStartDirZ", data->showerStartDirZ.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    
    // Angular info
    createArrayBranch(event_tree_, "showerCosTheta", data->showerCosTheta.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerCosThetaY", data->showerCosThetaY.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerDistToVtx", data->showerDistToVtx.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    
    // CNN scores
    createArrayBranch(event_tree_, "showerElScore", data->showerElScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerPhScore", data->showerPhScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerMuScore", data->showerMuScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerPiScore", data->showerPiScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerPrScore", data->showerPrScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerPID", data->showerPID.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "I");
    
    // CNN quality
    createArrayBranch(event_tree_, "showerComp", data->showerComp.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerPurity", data->showerPurity.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerProcess", data->showerProcess.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "I");
    
    // Origin scores
    createArrayBranch(event_tree_, "showerPrimaryScore", data->showerPrimaryScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerFromNeutralScore", data->showerFromNeutralScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerFromChargedScore", data->showerFromChargedScore.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    
    // Energy
    createArrayBranch(event_tree_, "showerRecoE", data->showerRecoE.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerRecoEU", data->showerRecoEU.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerRecoEV", data->showerRecoEV.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    createArrayBranch(event_tree_, "showerRecoEY", data->showerRecoEY.data(), 
                     "nShowers", EventData::MAX_SHOWERS, "F");
    
    // Truth matching (MC only)
    if (is_mc_) {
        createArrayBranch(event_tree_, "showerTruePDG", data->showerTruePDG.data(), 
                         "nShowers", EventData::MAX_SHOWERS, "I");
        createArrayBranch(event_tree_, "showerTrueTID", data->showerTrueTID.data(), 
                         "nShowers", EventData::MAX_SHOWERS, "I");
        createArrayBranch(event_tree_, "showerTrueMID", data->showerTrueMID.data(), 
                         "nShowers", EventData::MAX_SHOWERS, "I");
        createArrayBranch(event_tree_, "showerTrueE", data->showerTrueE.data(), 
                         "nShowers", EventData::MAX_SHOWERS, "F");
        createArrayBranch(event_tree_, "showerTrueComp", data->showerTrueComp.data(), 
                         "nShowers", EventData::MAX_SHOWERS, "F");
        createArrayBranch(event_tree_, "showerTruePurity", data->showerTruePurity.data(), 
                         "nShowers", EventData::MAX_SHOWERS, "F");
    }
}

void BranchManager::setupEventFeatureBranches(EventData* data) {
    // Total reconstructed energy
    event_tree_->Branch("recoNuE", &data->recoNuE, "recoNuE/F");
    
    // PCA analysis
    event_tree_->Branch("eventPCAxis0", data->eventPCAxis0.data(), "eventPCAxis0[3]/F");
    event_tree_->Branch("eventPCAxis1", data->eventPCAxis1.data(), "eventPCAxis1[3]/F");
    event_tree_->Branch("eventPCAxis2", data->eventPCAxis2.data(), "eventPCAxis2[3]/F");
    event_tree_->Branch("eventPCEigenVals", data->eventPCEigenVals.data(), "eventPCEigenVals[3]/F");
    event_tree_->Branch("eventPCProjMaxGap", &data->eventPCProjMaxGap, "eventPCProjMaxGap/F");
    event_tree_->Branch("eventPCProjMaxDist", &data->eventPCProjMaxDist, "eventPCProjMaxDist/F");
}

template<typename T>
void BranchManager::createArrayBranch(TTree* tree, const char* name, T* address, 
                                     const char* size_var, int max_size, const char* type) {
    std::stringstream ss;
    ss << name << "[" << size_var << "]/" << type;
    tree->Branch(name, address, ss.str().c_str());
}

void BranchManager::fillEventTree() {
    if (event_tree_) {
        event_tree_->Fill();
    }
}

void BranchManager::fillPOTTree() {
    if (pot_tree_) {
        pot_tree_->Fill();
    }
}

void BranchManager::write() {
    if (output_file_) {
        output_file_->cd();
        
        if (event_tree_) {
            event_tree_->Write();
        }
        
        if (pot_tree_) {
            pot_tree_->Write();
        }
    }
}

Long64_t BranchManager::getEventEntries() const {
    if (event_tree_) {
        return event_tree_->GetEntries();
    }
    return 0;
}

Long64_t BranchManager::getPOTEntries() const {
    if (pot_tree_) {
        return pot_tree_->GetEntries();
    }
    return 0;
}

} // namespace gen2ntuple