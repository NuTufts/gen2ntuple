#pragma once

#include "EventData.h"
#include <vector>
#include <memory>

// Forward declarations
namespace larlite {
    class storage_manager;
    class shower;
    class pcaxis;
}

namespace larcv {
    class IOManager;
}

namespace gen2ntuple {

class ShowerProcessor {
public:
    ShowerProcessor();
    ~ShowerProcessor() = default;
    
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
    
    // Shower processing methods
    bool extractShowerInfo(larlite::storage_manager* larlite_io,
                          larcv::IOManager* larcv_io,
                          EventData* event_data);
    
    bool calculateShowerGeometry(const larlite::shower& shower, int shower_idx,
                                EventData* event_data);
    
    bool calculateShowerAngles(const larlite::shower& shower, int shower_idx,
                              EventData* event_data);
    
    bool calculateShowerCharge(larlite::storage_manager* larlite_io,
                              larcv::IOManager* larcv_io,
                              int shower_idx, EventData* event_data);
    
    bool calculateShowerEnergy(const larlite::shower& shower, int shower_idx,
                              EventData* event_data);
    
    // Truth matching (MC only)
    bool performTruthMatching(larlite::storage_manager* larlite_io,
                             int shower_idx, EventData* event_data);
    
    // Geometric calculations
    float calculateDistanceToVertex(const larlite::shower& shower) const;
    float calculateCosTheta(const larlite::shower& shower) const;
    float calculateCosThetaY(const larlite::shower& shower) const;
    std::vector<float> getShowerDirection(const larlite::shower& shower) const;
    
    // Energy reconstruction
    float calculateElectronEnergy(const larlite::shower& shower) const;
    float calculatePhotonEnergy(const larlite::shower& shower) const;
    float calculateShowerEnergyByPlane(const larlite::shower& shower, int plane) const;
    
    // Secondary shower identification
    bool isSecondaryShower(const larlite::shower& shower) const;
    
    // Constants
    static constexpr float BEAM_DIR_X = 0.0f;
    static constexpr float BEAM_DIR_Y = 0.0f; 
    static constexpr float BEAM_DIR_Z = 1.0f;
    
    static constexpr float GRAVITY_DIR_X = 0.0f;
    static constexpr float GRAVITY_DIR_Y = -1.0f;
    static constexpr float GRAVITY_DIR_Z = 0.0f;
    
    // Particle ID constants for EM showers
    static constexpr int PID_UNKNOWN = -1;
    static constexpr int PID_ELECTRON = 11;
    static constexpr int PID_PHOTON = 22;
    static constexpr int PID_PI0 = 111;
    
    // Energy constants
    static constexpr float MEV_TO_GEV = 0.001f;
    static constexpr float ELECTRON_MASS = 0.511f; // MeV
};

} // namespace gen2ntuple
