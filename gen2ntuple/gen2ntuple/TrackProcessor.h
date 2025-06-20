#pragma once

#include "EventData.h"
#include "RecoData.h"
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

namespace larpid {
namespace model {
    class TorchModel;
}
}

#include "larpid/data/CropPixData_t.h"
#include "ublarcvapp/MCTools/MCPixelPGraph.h"
#include "ublarcvapp/MCTools/MCPixelPMap.h"

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
    // Set ProngCNN model (non-owning pointer - caller retains ownership)
    void setProngCNNInterface(larpid::model::TorchModel* larpid_model ) {
        larpid_cnn_ = larpid_model;
    }
    
    void setMCPixelUtils( ublarcvapp::mctools::MCPixelPGraph* mcpg,
                          ublarcvapp::mctools::MCPixelPMap* mcpm )
    {
        _mcpg = mcpg;
        _mcpm = mcpm;
    }

    // Main processing method
    bool processEvent(larlite::storage_manager* larlite_io,
                     larcv::IOManager* larcv_io,
                     EventData* event_data,
                     RecoData* reco_data);
    
private:
    bool is_mc_;
    float vertex_x_, vertex_y_, vertex_z_;
    larpid::model::TorchModel* larpid_cnn_;  // Non-owning pointer
    ublarcvapp::mctools::MCPixelPGraph* _mcpg;
    ublarcvapp::mctools::MCPixelPMap*   _mcpm;
    
    // Track processing methods
    bool extractTrackInfo( larlite::storage_manager* larlite_io,
                           larcv::IOManager* larcv_io,
                           EventData* event_data,
                           RecoData* reco_data );
    
    bool calculateTrackGeometry(const larlite::track& track, 
                                const larlite::larflowcluster& cluster,
                                int track_idx,
                                EventData* event_data);
    
    bool calculateTrackAngles(const larlite::track& track, int track_idx,
                             EventData* event_data);
    
    bool calculateTrackCharge( larcv::IOManager* larcv_io,
                               EventData* event_data,
                               RecoData* reco_data,
                               int vtxIdx,
                               int track_idx );
    
    bool calculateTrackEnergy(const larlite::track& track, int track_idx,
                             EventData* event_data);
    
    // Truth matching (MC only)
    bool performTruthMatching(larlite::storage_manager* larlite_io,
                             int track_idx, EventData* event_data);
    
    // Geometric calculations
    float calculateTrackLength(const larlite::track& track) const;
    float calculateDistanceToVertex(const larlite::track& track) const;
    std::vector<float> getTrackDirection(const larlite::track& track) const;
    
    // Energy reconstruction
    float calculateMuonEnergy(const larlite::track& track) const;
    float calculatePionRangeEnergy(float track_length) const;
    
    // ProngCNN integration
    bool runProngCNN(larcv::IOManager* larcv_io,
                     const larlite::track& track,
                     const larlite::larflowcluster& cluster,
                     int track_idx,
                     int& num_good_planes,
                     EventData* event_data,
                     RecoData* reco_data);
    
    void setDefaultPIDScores(int track_idx, EventData* event_data);
    
    bool updateEnergyBasedOnPID( int vtxIdx, int track_idx,
                                 EventData* event_data,
                                 RecoData* reco_data );
    
    float calculateProtonEnergy(const larlite::track& track) const;
    
    int getPIDFromScores(float el_score, float ph_score, float mu_score,
                        float pi_score, float pr_score) const;

    bool getMCProngParticles( larcv::IOManager* larcv_io,
            std::vector< std::vector<larpid::data::CropPixData_t> >& prong_vv,
            EventData* event_data, int track_idx );
    
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
