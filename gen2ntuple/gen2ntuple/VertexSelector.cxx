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
                                 std::vector<larflow::reco::NuVertexCandidate>* nuvtx_v,
                                 EventData* event_data) {
    
    if (!event_data) {
        std::cerr << "VertexSelector: EventData pointer is null" << std::endl;
        return false;
    }

    // Calculate Flash Prediction to each vertex
    calculateFlashPredictions(larlite_io, larcv_io, nuvtx_v);
    
    // Find the best vertex candidate
    if (!findBestVertex(larlite_io, larcv_io, event_data)) {
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
    if (!calculateVertexQuality(larlite_io, larcv_io, event_data)) {
        return false;
    }
    
    // Process keypoints if enabled
    if (include_keypoints_) {
        if (!processKeypoints(larcv_io, event_data)) {
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
                                   EventData* event_data) {
    
    // Get vertex collection from LArLite
    auto ev_vertex = larlite_io->get_data<larlite::event_vertex>("nuvertex");
    if (!ev_vertex || ev_vertex->size() == 0) {
        return false;
    }
    
    // For now, select the first vertex (in practice, this would be more sophisticated)
    // The Python code uses various selection criteria based on scores, cosmic rejection, etc.
    const auto& vertex = ev_vertex->at(0);
    
    event_data->foundVertex = true;
    event_data->vtxX = vertex.X();
    event_data->vtxY = vertex.Y(); 
    event_data->vtxZ = vertex.Z();
    
    // Get vertex score if available (would need to be stored in vertex object)
    event_data->vtxScore = 1.0f; // Placeholder
    
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
                                           EventData* event_data) {
    
    // Placeholder values - in practice these would be calculated from:
    // - Pixel analysis around vertex
    // - Cosmic track association  
    // - Keypoint analysis
    
    event_data->vtxKPtype = 0;
    event_data->vtxKPscore = 0.0f;
    event_data->vtxMaxIntimePixelSum = 1000.0f;
    event_data->vtxFracHitsOnCosmic = 0.1f;
    
    // Pixel fractions (placeholder)
    for (int i = 0; i < 3; i++) {
        event_data->fracUnrecoIntimePixels[i] = 0.05f;
        event_data->fracRecoOuttimePixels[i] = 0.02f;
    }
    
    return true;
}

bool VertexSelector::processKeypoints(larcv::IOManager* larcv_io, EventData* event_data) {
    
    // Get keypoint data from LArCV
    // This is a simplified version - the full implementation would process
    // the keypoint pixel clusters and extract features
    
    event_data->nKeypoints = 0; // Placeholder
    
    // In practice, would iterate through keypoint clusters:
    // - Extract cluster properties
    // - Calculate scores and positions
    // - Fill arrays up to MAX_KEYPOINTS
    
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
                                               std::vector<larflow::reco::NuVertexCandidate>* nuvtx_v) 
{

    _nuvtx_sinkhorn_div_v.clear();
    _nuvtx_flashpred_v.clear();

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
        }       
    }
    else {
        // Create flat dummy opflash
        obs_total_pe = 0.0;
        for (int pmt = 0; pmt < 32; pmt++) {
            obs_pe_per_pmt[pmt] = 1.0/32.0;
        }        
    }

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

} // namespace gen2ntuple
