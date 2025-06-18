#pragma once

#include "EventData.h"
#include <vector>
#include <memory>
#include <map>
#include <tuple>

#include "TVector3.h"

// Forward declarations
namespace larlite {
    class storage_manager;
    class mctruth;
    class mcpart;
}
namespace larutil {
    class SpaceChargeMicroBooNE;
}
namespace larcv {
    class IOManager;
}
namespace ublarcvapp {
namespace mctools {
    class MCPixelPGraph;
}
}


namespace gen2ntuple {

class MCTruthProcessor {
public:
    MCTruthProcessor();
    ~MCTruthProcessor();
    
    // Main processing method
    bool processEvent(larlite::storage_manager* larlite_io, 
                      larcv::IOManager* larcv_io,
                      EventData* event_data);
    
    // Cross-section weight handling
    void setWeightFile(const std::string& weight_file) { weight_file_ = weight_file; }
    bool loadWeights();
    float getEventWeight(int run, int subrun, int event) const;
    
private:
    std::string weight_file_;
    
    // Processing methods
    bool extractNeutrinoTruth(larlite::storage_manager* larlite_io, EventData* event_data);
    bool extractTrueVertex(larlite::storage_manager* larlite_io, EventData* event_data);
    bool extractPrimaryParticles(larlite::storage_manager* larlite_io, EventData* event_data);
    bool extractSimulatedParticles(larlite::storage_manager* larlite_io, EventData* event_data);
    
    // Helper methods
    bool isNeutrinoInteraction(const larlite::mctruth& mctruth) const;
    bool isPrimaryParticle(const larlite::mcpart& mcpart) const;
    int  getParticleProcess(const larlite::mcpart& mcpart) const;
    TVector3 getSCECorrectedPos( double x, double y, double z );
    
    // Weight storage - uses pickle loading via Python C-API
    std::map<std::tuple<int,int,int>, float> event_weights_;

    // SpaceCharge Utility
    larutil::SpaceChargeMicroBooNE* pSCE;

    // Particle Graph Utility
    ublarcvapp::mctools::MCPixelPGraph* pmcpg;
    
    // Python C-API handling
    bool initializePython();
    void finalizePython();
    bool python_initialized_ = false;
    
    // Constants
    static constexpr int NU_MUON_PDG = 14;
    static constexpr int NU_ELECTRON_PDG = 12;
    static constexpr int MUON_PDG = 13;
    static constexpr int ELECTRON_PDG = 11;
    static constexpr int PHOTON_PDG = 22;
    static constexpr int PROTON_PDG = 2212;
    static constexpr int NEUTRON_PDG = 2112;
    static constexpr int PION_PLUS_PDG = 211;
    static constexpr int PION_MINUS_PDG = -211;
    static constexpr int PION_ZERO_PDG = 111;
};

} // namespace gen2ntuple