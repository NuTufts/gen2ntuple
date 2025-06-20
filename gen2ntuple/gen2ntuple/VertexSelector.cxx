#include "VertexSelector.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// LArLite includes
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/vertex.h"
#include "larlite/DataFormat/mctruth.h"
#include "larlite/DataFormat/opflash.h"
#include "larlite/LArUtil/SpaceChargeMicroBooNE.h"

// LArCV includes
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventPixel2D.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// UBLArCVApp includes (for truth matching)
#include "ublarcvapp/MCTools/NeutrinoVertex.h"

#include "larflow/Reco/NuVertexFlashPrediction.h"
#include "larflow/Reco/SinkhornFlashDivergence.h"

// Include concrete vertex selection classes
#include "SelectNuVertex.h"
#include "HighestKPRankWithVisEnergy.h"
#include "GetNuCandidateIntimeCharge.h"
#include "FlashPEVertexSelection.h"

namespace gen2ntuple {

VertexSelector::VertexSelector() 
    : is_mc_(false), include_keypoints_(false) 
{
    _nuvtx_flashpred_v.clear();

    predictor = new larflow::reco::NuVertexFlashPrediction;
    // Configure flash predictor with standard parameters
    predictor->setChargeToPhotonParams(
        200.0,    // adc_per_electron
        23.6e-3,  // mev_per_electron (MeV)
        24000.0,  // photons_per_mev
        0.7       // recombination_factor
    );
    predictor->setTrackConversionParams(
        3,      // dcol
        3,      // drow
        0.3,    // minstepsize (cm)
        0.5     // maxstepsize (cm)
    );
    predictor->setShowerConversionParams(
        3,      // dcol
        3       // drow
    );     

    sinkhorn_calc = new larflow::reco::SinkhornFlashDivergence;
    
    // Load all vertex selection methods
    loadVertexSelectionMethods();
}

VertexSelector::~VertexSelector()
{
    // Clear flash prediction vector
    delete predictor;
    predictor = nullptr;
    delete sinkhorn_calc;
    sinkhorn_calc = nullptr;
    _nuvtx_flashpred_v.clear();
}

bool VertexSelector::processEvent(larlite::storage_manager* larlite_io,
                                 larcv::IOManager* larcv_io,
                                 EventData* event_data,
                                 RecoData* reco_data,
                                 std::string vertex_selector ) {
    
    if (!event_data) {
        std::cerr << "VertexSelector: EventData pointer is null" << std::endl;
        return false;
    }
    
    if (!reco_data || !reco_data->nuvtx_v) {
        std::cerr << "VertexSelector: RecoData or vertex candidates not available" << std::endl;
        return false;
    }

    // Calculate Flash Prediction to each vertex
    calculateFlashPredictions(larlite_io, larcv_io, event_data, reco_data);
    
    // Find the best vertex candidate
    if (!findBestVertex(larlite_io, larcv_io, event_data, reco_data, vertex_selector )) {
        // No vertex found - set defaults
        event_data->foundVertex = false;
        event_data->vtxScore = -1.0f;
        event_data->vtxX = event_data->vtxY = event_data->vtxZ = -999.0f;
        event_data->vtxIsFiducial = false;
        event_data->vtxContainment = 0.0f;
        
        if (is_mc_) {
            event_data->vtxDistToTrue = -1.0f;
        }
        
        return true; // Not finding a vertex is OK
    }
    
    // Extract vertex information
    if (!extractVertexInfo(larlite_io, event_data)) {
        return false;
    }
    
    // Calculate vertex quality metrics
    if (!calculateVertexQuality( larlite_io, larcv_io, 
                                 event_data, reco_data )) {
        return false;
    }
    
    // Process keypoints if enabled
    if (include_keypoints_) {
        if (!processKeypoints(reco_data, event_data)) {
            return false;
        }
    }
    
    // Calculate truth distance for MC
    if (is_mc_) {
        if (!calculateTruthDistance(larlite_io, event_data)) {
            return false;
        }
    }
    
    return true;
}

bool VertexSelector::findBestVertex(larlite::storage_manager* larlite_io,
                                   larcv::IOManager* larcv_io,
                                   EventData* event_data,
                                   RecoData* reco_data,
                                   std::string vertex_selector ) {
    
    // Check if we have vertex candidates
    if (!reco_data || !reco_data->nuvtx_v || reco_data->nuvtx_v->empty()) {
        std::cout << "VertexSelector: No vertex candidates available" << std::endl;
        event_data->foundVertex = false;
        return false;
    }
    
    // Find the requested selection method
    auto method_iter = _vertex_selection_methods_m.find(vertex_selector);
    if (method_iter == _vertex_selection_methods_m.end()) {
        std::cerr << "VertexSelector: Unknown vertex selection method: " << vertex_selector << std::endl;
        std::cerr << "Available methods: ";
        for (const auto& pair : _vertex_selection_methods_m) {
            std::cerr << pair.first << " ";
        }
        std::cerr << std::endl;
        return false;
    }
    
    // Pass flash predictions to the selection method if it supports it
    method_iter->second->setFlashPredictions(_nuvtx_flashpred_v, _nuvtx_sinkhorn_div_v);
    
    // Use the selected method to find the best vertex
    int selected_index = method_iter->second->selectVertex(larcv_io, larlite_io, event_data, reco_data);
    
    if (selected_index < 0 || selected_index >= static_cast<int>(reco_data->nuvtx_v->size())) {
        std::cout << "VertexSelector: No vertex selected by method " << vertex_selector << std::endl;
        event_data->foundVertex = false;
        return false;
    }
    
    // Get the selected vertex
    const auto& selected_vtx = reco_data->nuvtx_v->at(selected_index);
    
    // Fill event data with selected vertex information
    event_data->foundVertex = true;
    event_data->vtxIndex = selected_index;
    event_data->vtxX = selected_vtx.pos[0];
    event_data->vtxY = selected_vtx.pos[1]; 
    event_data->vtxZ = selected_vtx.pos[2];
    event_data->vtxScore   = selected_vtx.netScore;
    event_data->vtxKPtype  = selected_vtx.keypoint_type;
    event_data->vtxKPscore = selected_vtx.netScore;

    if ( selected_index>=0 && selected_index<(int)_nuvtx_sinkhorn_div_v.size() ) {
        /// fill flash prediction variables
        event_data->predictedPEtotal = 0.;
        for (int ipmt=0; ipmt<32; ipmt++) {
            event_data->predictedPE[ipmt] = _nuvtx_flashpred_v[selected_index][ipmt];
            event_data->predictedPEtotal += _nuvtx_flashpred_v[selected_index][ipmt];
        }
        event_data->sinkhorn_div = _nuvtx_sinkhorn_div_v.at(selected_index);
        event_data->fracerrPE = (event_data->predictedPEtotal-event_data->observedPEtotal)/(0.1+event_data->observedPEtotal);
    }
    
    std::cout << "VertexSelector: Selected vertex " << selected_index 
              << " at (" << event_data->vtxX << ", " << event_data->vtxY 
              << ", " << event_data->vtxZ << ") with score " << event_data->vtxScore << std::endl;
    
    return true;
}

bool VertexSelector::extractVertexInfo(larlite::storage_manager* larlite_io,
                                      EventData* event_data) {
    
    // Check if vertex is in fiducial volume
    event_data->vtxIsFiducial = checkFiducialVolume(event_data->vtxX, 
                                                    event_data->vtxY, 
                                                    event_data->vtxZ);
    
    // Calculate containment
    event_data->vtxContainment = calculateContainment(event_data->vtxX,
                                                     event_data->vtxY,
                                                     event_data->vtxZ);
    
    return true;
}

bool VertexSelector::calculateVertexQuality(larlite::storage_manager* larlite_io,
                                           larcv::IOManager* larcv_io,
                                           EventData* event_data,
                                           RecoData* reco_data ) {
    
    int vtxIdx = event_data->vtxIndex;

    if ( vtxIdx<0 )
        return true;

    auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
    auto const& nusel = reco_data->nusel_v->at(vtxIdx);
    
    event_data->vtxFracHitsOnCosmic = nusel.frac_allhits_on_cosmic;
    
    if ( nusel.unreco_fraction_v.size()>=3 ) {
        for (int i = 0; i < 3; i++) {
            event_data->fracUnrecoIntimePixels[i] = nusel.unreco_fraction_v.at(i);
        }
    }

    if ( nusel.unreco_fraction_v.size()>=6 ) {
        for (int i=0; i<3; i++) {
            event_data->fracRecoOuttimePixels[i]  = nusel.unreco_fraction_v.at(3+i);
        }
    }
    
    return true;
}

bool VertexSelector::processKeypoints(RecoData* reco_data, 
                                      EventData* event_data) 
{     
    if ( reco_data->kpc_nu_v==nullptr ) {
        return true;
    }

    event_data->nKeypoints = 0;

    int iKP=0;

    std::vector< std::vector< larflow::reco::KPCluster >* > kpcontainers;
    kpcontainers.push_back( reco_data->kpc_nu_v );
    kpcontainers.push_back( reco_data->kpc_track_v );
    kpcontainers.push_back( reco_data->kpc_shower_v );
    kpcontainers.push_back( reco_data->kpc_cosmic_v );
    
    for ( auto&  pkp_v : kpcontainers ) {

        if ( iKP>= EventData::MAX_KEYPOINTS )
            break;

        for (auto const& kpc : *pkp_v ) {
            if ( iKP>=EventData::MAX_KEYPOINTS )
                break;

            event_data->kpClusterType[iKP] = kpc._cluster_type;
            event_data->kpFilterType[iKP]  = 0;
            event_data->kpMaxScore[iKP]    = kpc.max_score;
            event_data->kpMaxPosX[iKP]     = kpc.max_pt_v[0];
            event_data->kpMaxPosY[iKP]     = kpc.max_pt_v[1];
            event_data->kpMaxPosZ[iKP]     = kpc.max_pt_v[2];
            iKP++;
        }
    }
    event_data->nKeypoints = iKP;
    
    return true;
}

bool VertexSelector::checkFiducialVolume(float x, float y, float z) const {
    
    // Wire Cell fiducial volume definition: >3cm from SCE-corrected edges
    bool x_fid = (x > DETECTOR_X_MIN + FV_BORDER) && (x < DETECTOR_X_MAX - FV_BORDER);
    bool y_fid = (y > DETECTOR_Y_MIN + FV_BORDER) && (y < DETECTOR_Y_MAX - FV_BORDER);  
    bool z_fid = (z > DETECTOR_Z_MIN + FV_BORDER) && (z < DETECTOR_Z_MAX - FV_BORDER);
    
    return x_fid && y_fid && z_fid;
}

float VertexSelector::calculateContainment(float x, float y, float z) const {
    
    // Calculate minimum distance to detector boundary
    float dx_min = std::min(x - DETECTOR_X_MIN, DETECTOR_X_MAX - x);
    float dy_min = std::min(y - DETECTOR_Y_MIN, DETECTOR_Y_MAX - y);
    float dz_min = std::min(z - DETECTOR_Z_MIN, DETECTOR_Z_MAX - z);
    
    float min_dist = std::min({dx_min, dy_min, dz_min});
    
    // Convert to containment score (larger distance = better containment)
    return std::max(0.0f, min_dist);
}

bool VertexSelector::calculateTruthDistance(larlite::storage_manager* larlite_io,
                                           EventData* event_data) {
    
    if (!is_mc_) {
        return true;
    }
    
    // Get MC truth information
    auto ev_mctruth = larlite_io->get_data<larlite::event_mctruth>("generator");
    if (!ev_mctruth || ev_mctruth->size() == 0) {
        event_data->vtxDistToTrue = -1.0f;
        return true;
    }
    
    try {
        // Use NeutrinoVertex to get SCE-corrected position
        auto mc_nu_vertexer = std::make_unique<ublarcvapp::mctools::NeutrinoVertex>();
        auto sce = std::make_unique<larutil::SpaceChargeMicroBooNE>();
        
        std::vector<float> mc_nu_vertex = mc_nu_vertexer->getPos3DwSCE(*larlite_io, sce.get());
        
        if (mc_nu_vertex.size() >= 3) {
            float dx = event_data->vtxX - mc_nu_vertex[0];
            float dy = event_data->vtxY - mc_nu_vertex[1]; 
            float dz = event_data->vtxZ - mc_nu_vertex[2];
            
            event_data->vtxDistToTrue = std::sqrt(dx*dx + dy*dy + dz*dz);
        } else {
            event_data->vtxDistToTrue = -1.0f;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "VertexSelector: Error calculating truth distance: " << e.what() << std::endl;
        event_data->vtxDistToTrue = -1.0f;
    }
    
    return true;
}

int VertexSelector::calculateFlashPredictions( larlite::storage_manager* larlite_io,
                                               larcv::IOManager* larcv_io,
                                               EventData* event_data,
                                               RecoData* reco_data )

{

    _nuvtx_sinkhorn_div_v.clear();
    _nuvtx_flashpred_v.clear();

    std::vector<larflow::reco::NuVertexCandidate>* nuvtx_v = reco_data->nuvtx_v;

    // Get Wire Plane Images
    auto ev_img = (larcv::EventImage2D*)(larcv_io->get_data(larcv::kProductImage2D, "wire"));
    if (!ev_img || ev_img->Image2DArray().size() < 3) {
        throw std::runtime_error("VertexSelector::calculateFlashPredictions - Wrong number of images");
    }
    auto const& adc_v = ev_img->as_vector();
    float adc_threshold = 10.0;

    // Get Observed Flash
    auto ev_opflash = (larlite::event_opflash*)(larlite_io->get_data(larlite::data::kOpFlash, "simpleFlashBeam"));    
    bool has_flash = (ev_opflash && ev_opflash->size() > 0);  
    std::vector<float> obs_pe_per_pmt(32,0.0);
    float obs_total_pe = 0.;
    if ( has_flash ) {
        const auto& flash = ev_opflash->at(0);
        obs_total_pe = flash.TotalPE();
        for (int pmt = 0; pmt < 32; pmt++) {
            obs_pe_per_pmt[pmt] = flash.PE(pmt);
            event_data->observedPE[pmt] = flash.PE(pmt);
        }       
    }
    else {
        // Create flat dummy opflash
        obs_total_pe = 0.0;
        for (int pmt = 0; pmt < 32; pmt++) {
            obs_pe_per_pmt[pmt] = 1.0/32.0;
            event_data->observedPE[pmt] = 1.0/32.0;
            obs_total_pe += 1.0/32.0;
        }        
    }
    
    // store into EventData
    event_data->observedPEtotal = obs_total_pe;

    // calculate prediction vs. observe metrics

    float sinkhorn_regularization = 1.0; 

    for (auto& vtx : *nuvtx_v ) {

        // Make Prediction Flash
        auto predicted_flash_all = predictor->predictFlash(
            vtx,
            adc_v,
            adc_threshold,
            true,   // use_trilinear
            false   // primary_prongs_only = false (all particles)
        );

        std::vector<float> pred_pe_per_pmt(32,0.0);
        for (int pmt = 0; pmt < 32; pmt++) {
            pred_pe_per_pmt[pmt] = predicted_flash_all.PE(pmt);
        }        

        _nuvtx_flashpred_v.push_back( pred_pe_per_pmt );

        // Calculate Sinkhorn Divergence predicted and observed
        float sinkhorn_div = sinkhorn_calc->calculateDivergence(
            pred_pe_per_pmt,
            obs_pe_per_pmt,
            sinkhorn_regularization,
            100,
            1e-6
        );
        _nuvtx_sinkhorn_div_v.push_back( sinkhorn_div );
    }

    return 0;
}

bool VertexSelector::loadVertexSelectionMethods() {
    // Clear any existing methods
    _vertex_selection_methods_m.clear();
    
    // Register all available vertex selection methods
    _vertex_selection_methods_m["select_nu_vertex"] = 
        std::make_unique<SelectNuVertex>();
        
    _vertex_selection_methods_m["highest_kprank_with_visenergy"] = 
        std::make_unique<HighestKPRankWithVisEnergy>();
        
    _vertex_selection_methods_m["highest_intime_reco_frac"] = 
        std::make_unique<GetNuCandidateIntimeCharge>();
        
    _vertex_selection_methods_m["flash_pe_selection"] = 
        std::make_unique<FlashPEVertexSelection>();
    
    std::cout << "VertexSelector: Loaded " << _vertex_selection_methods_m.size() 
              << " vertex selection methods" << std::endl;
    
    return true;
}

} // namespace gen2ntuple
