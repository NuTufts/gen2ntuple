#include "VertexSelector.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// LArLite includes
#include "DataFormat/storage_manager.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"

// LArCV includes
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventPixel2D.h"

// UBLArCVApp includes (for truth matching)
#include "ublarcvapp/MCTools/NeutrinoVertex.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"

namespace gen2ntuple {

VertexSelector::VertexSelector() 
    : is_mc_(false), include_keypoints_(false) {
}

bool VertexSelector::processEvent(larlite::storage_manager* larlite_io,
                                 larcv::IOManager* larcv_io,
                                 EventData* event_data) {
    
    if (!event_data) {
        std::cerr << "VertexSelector: EventData pointer is null" << std::endl;
        return false;
    }
    
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

} // namespace gen2ntuple
