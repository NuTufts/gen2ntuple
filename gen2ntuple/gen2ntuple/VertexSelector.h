#pragma once

#include "EventData.h"
#include <vector>
#include <memory>

// Forward declarations
namespace larlite {
    class storage_manager;
}

namespace larcv {
    class IOManager;
}

namespace gen2ntuple {

class VertexSelector {
public:
    VertexSelector();
    ~VertexSelector() = default;
    
    // Configuration
    void setMCMode(bool is_mc) { is_mc_ = is_mc; }
    void setKeypointMode(bool include_keypoints) { include_keypoints_ = include_keypoints; }
    
    // Main processing method
    bool processEvent(larlite::storage_manager* larlite_io,
                     larcv::IOManager* larcv_io,
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
};

} // namespace gen2ntuple
