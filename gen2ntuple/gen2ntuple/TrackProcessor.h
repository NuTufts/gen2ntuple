#pragma once

#include "EventData.h"
#include <vector>
#include <memory>

// Forward declarations
namespace larlite {
    class storage_manager;
    class track;
    class pcaxis;
}

namespace larcv {
    class IOManager;
}

namespace gen2ntuple {

class TrackProcessor {
public:
    TrackProcessor();
    ~TrackProcessor() = default;
    
    // Configuration
    void setMCMode(bool is_mc) { is_mc_ = is_mc; }
    void setVertexPosition(float x, float y, float z) { 
        vertex_x_ = x; vertex_y_ = y; vertex_z_ = z; 
    }
    
    // Main processing method
    bool processEvent(larlite::storage_manager* larlite_io,
                     larcv::IOManager* larcv_io,
                     EventData* event_data);
    
private:
    bool is_mc_;
    float vertex_x_, vertex_y_, vertex_z_;
    
    // Track processing methods
    bool extractTrackInfo(larlite::storage_manager* larlite_io,
                         EventData* event_data);
    
    bool calculateTrackGeometry(const larlite::track& track, int track_idx,
                               EventData* event_data);
    
    bool calculateTrackAngles(const larlite::track& track, int track_idx,
                             EventData* event_data);
    
    bool calculateTrackCharge(larlite::storage_manager* larlite_io,
                             larcv::IOManager* larcv_io,
                             int track_idx, EventData* event_data);
    
    bool calculateTrackEnergy(const larlite::track& track, int track_idx,
                             EventData* event_data);
    
    // Truth matching (MC only)
    bool performTruthMatching(larlite::storage_manager* larlite_io,
                             int track_idx, EventData* event_data);
    
    // Geometric calculations
    float calculateTrackLength(const larlite::track& track) const;
    float calculateDistanceToVertex(const larlite::track& track) const;
    float calculateCosTheta(const larlite::track& track) const;
    float calculateCosThetaY(const larlite::track& track) const;
    std::vector<float> getTrackDirection(const larlite::track& track) const;
    
    // Energy reconstruction
    float calculateMuonEnergy(const larlite::track& track) const;
    float calculatePionRangeEnergy(float track_length) const;
    
    // Secondary track identification
    bool isSecondaryTrack(const larlite::track& track) const;
    
    // Constants
    static constexpr float BEAM_DIR_X = 0.0f;
    static constexpr float BEAM_DIR_Y = 0.0f; 
    static constexpr float BEAM_DIR_Z = 1.0f;
    
    static constexpr float GRAVITY_DIR_X = 0.0f;
    static constexpr float GRAVITY_DIR_Y = -1.0f;
    static constexpr float GRAVITY_DIR_Z = 0.0f;
    
    // Particle ID constants
    static constexpr int PID_UNKNOWN = -1;
    static constexpr int PID_ELECTRON = 11;
    static constexpr int PID_MUON = 13;
    static constexpr int PID_PION = 211;
    static constexpr int PID_PROTON = 2212;
    static constexpr int PID_PHOTON = 22;
};

} // namespace gen2ntuple
