#pragma once

#include "EventData.h"
#include <vector>
#include <memory>

#include "larflow/Reco/NuVertexCandidate.h"

// Forward declarations
namespace larlite {
    class storage_manager;
}

namespace larcv {
    class IOManager;
}

namespace larflow {
namespace reco {
    class NuVertexFlashPrediction;
    class SinkhornFlashDivergence;
}
}

namespace gen2ntuple {

class VertexSelector {
public:
    VertexSelector();
    ~VertexSelector();
    
    // Configuration
    void setMCMode(bool is_mc) { is_mc_ = is_mc; }
    void setKeypointMode(bool include_keypoints) { include_keypoints_ = include_keypoints; }
    
    // Main processing method
    bool processEvent(larlite::storage_manager* larlite_io,
                     larcv::IOManager* larcv_io,
                     std::vector<larflow::reco::NuVertexCandidate>* nuvtx_v,
                     EventData* event_data);
    
private:
    bool is_mc_;
    bool include_keypoints_;
    
    // Vertex selection methods
    bool findBestVertex(larlite::storage_manager* larlite_io,
                       larcv::IOManager* larcv_io,
                       EventData* event_data);
    
    bool extractVertexInfo(larlite::storage_manager* larlite_io,
                          EventData* event_data);
    
    bool calculateVertexQuality(larlite::storage_manager* larlite_io,
                               larcv::IOManager* larcv_io,
                               EventData* event_data);

    int calculateFlashPredictions( larlite::storage_manager* larlite_io,
                                   larcv::IOManager* larcv_io,
                                   std::vector<larflow::reco::NuVertexCandidate>* nuvtx_v);

    
    // Keypoint processing
    bool processKeypoints(larcv::IOManager* larcv_io, EventData* event_data);
    
    // Geometric calculations
    bool checkFiducialVolume(float x, float y, float z) const;
    float calculateContainment(float x, float y, float z) const;
    
    // Truth matching (MC only)
    bool calculateTruthDistance(larlite::storage_manager* larlite_io,
                                EventData* event_data);
    
    // Constants for fiducial volume (Wire Cell definition)
    static constexpr float FV_BORDER = 3.0; // cm
    static constexpr float DETECTOR_X_MIN = 0.0;
    static constexpr float DETECTOR_X_MAX = 256.35;
    static constexpr float DETECTOR_Y_MIN = -116.5;
    static constexpr float DETECTOR_Y_MAX = 116.5;
    static constexpr float DETECTOR_Z_MIN = 0.0;
    static constexpr float DETECTOR_Z_MAX = 1036.8;

    // Flash Prediction per vertex, for selection and also for EventData variables
    larflow::reco::NuVertexFlashPrediction* predictor;
    larflow::reco::SinkhornFlashDivergence* sinkhorn_calc;
    std::vector< std::vector<float> > _nuvtx_flashpred_v; ///< predictions for all nuvtx candidates
    std::vector< float > _nuvtx_sinkhorn_div_v;
};

} // namespace gen2ntuple
