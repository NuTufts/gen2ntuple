#pragma once

#include "EventData.h"
#include "RecoData.h"
#include <vector>
#include <memory>

// Forward declarations
namespace larlite {
    class storage_manager;
    class shower;
    class track;
    class pcaxis;
    class larflowcluster;
}

namespace larflow {
namespace reco {
    class NuVertexCandidate;
}
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

class ShowerProcessor {
public:
    ShowerProcessor();
    ~ShowerProcessor() = default;
    
    // Configuration
    void setMCMode(bool is_mc) { is_mc_ = is_mc; }
    void setVertexPosition(float x, float y, float z) { 
        vertex_x_ = x; vertex_y_ = y; vertex_z_ = z; 
    }
    // Set LArPID model (non-owning pointer - caller retains ownership)
    void setLArPIDInterface(larpid::model::TorchModel* larpid_model) {
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

    // Shower processing methods
    bool extractShowerInfo(larlite::storage_manager* larlite_io,
                          larcv::IOManager* larcv_io,
                          EventData* event_data,
                          RecoData* reco_data);
    
    bool calculateShowerGeometry(const larlite::track& shower_trunk, 
                                const larlite::larflowcluster& cluster,
                                int shower_idx,
                                EventData* event_data);
    
    bool calculateShowerAngles(const larlite::track& shower_trunk, int shower_idx,
                              EventData* event_data);
    
    bool calculateShowerCharge(larcv::IOManager* larcv_io,
                              EventData* event_data,
                              RecoData* reco_data,
                              int vtxIdx,
                              int shower_idx);
    
    bool calculateShowerEnergy(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx,
                              EventData* event_data);
    
    bool getMCProngParticles( larcv::IOManager* larcv_io,
            std::vector< std::vector<larpid::data::CropPixData_t> >& prong_vv,
            EventData* event_data, int shower_idx );
    
    // Geometric calculations
    float calculateDistanceToVertex(const larlite::track& shower_trunk) const;
    float calculateCosTheta(const larlite::track& shower_trunk) const;
    float calculateCosThetaY(const larlite::track& shower_trunk) const;
    std::vector<float> getShowerDirection(const larlite::track& shower_trunk) const;
    
    // Energy reconstruction
    float calculateElectronEnergy(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const;
    float calculatePhotonEnergy(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const;
    float calculateShowerEnergyByPlane(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx, int plane) const;
    
    // LArPID integration
    bool runLArPID(larcv::IOManager* larcv_io,
                   const larlite::track& shower_trunk,
                   const larlite::larflowcluster& cluster,
                   int shower_idx,
                   int& num_good_planes,
                   EventData* event_data,
                   RecoData* reco_data);
    
    void setDefaultPIDScores(int shower_idx, EventData* event_data);
    
    void updateEnergyBasedOnPID(const larflow::reco::NuVertexCandidate& nuvtx,
                                int shower_idx,
                                EventData* event_data);
    
    void calculateChargeFractions(larcv::IOManager* larcv_io,
                                  int shower_idx,
                                  EventData* event_data,
                                  RecoData* reco_data);
    
    int getPIDFromScores(float el_score, float ph_score, float mu_score,
                        float pi_score, float pr_score) const;
    
    // Secondary shower identification
    bool isSecondaryShower(const larlite::track& shower_trunk, const larlite::larflowcluster& cluster) const;
    
    // Additional shower calculation methods
    float calculateShowerLength(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const;
    float calculateOpeningAngle(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const;
    
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
    
    // Shower energy constants
    static constexpr float SHOWER_CALIB_FACTOR = 1.0f; // dE/dx calibration
    static constexpr float MIN_SHOWER_ENERGY = 10.0f;  // MeV minimum threshold
};

} // namespace gen2ntuple
