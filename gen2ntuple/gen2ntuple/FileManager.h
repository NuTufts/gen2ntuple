#pragma once

#include <string>
#include <vector>
#include <memory>

#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/NuSelectionVariables.h"

// Forward declarations
class TFile;
class TTree;
class TChain;

namespace larlite {
    class storage_manager;
}

namespace larcv {
    class IOManager;
}

namespace gen2ntuple {

class FileManager {
public:
    FileManager();
    ~FileManager();
    
    // Configuration
    bool setInputFiles(const std::vector<std::string>& kpsreco_files, 
                      const std::string& truth_file);
    void setMCMode(bool is_mc) { is_mc_ = is_mc; }
    void setDLAnaMode(bool is_dlana) { is_dlana_ = is_dlana; }
    
    // File operations
    bool openFiles();
    void closeFiles();
    bool nextEvent();
    void reset();
    
    // Access methods
    int getCurrentEntry() const { return current_entry_; }
    int getTotalEntries() const { return total_entries_; }
    bool isEndOfFile() const { return current_entry_ >= total_entries_; }
    
    // Data access
    larlite::storage_manager* getLarliteIO() { return larlite_io_.get(); }
    larcv::IOManager* getLarcvIO() { return larcv_io_.get(); }
    TChain* getRecoIO() { return kpsreco_.get(); }

    // Event synchronization
    struct EventID {
        int run;
        int subrun;
        int event;
        
        bool operator==(const EventID& other) const {
            return run == other.run && subrun == other.subrun && event == other.event;
        }
    };
    
    EventID getCurrentEventID() const;
    bool synchronizeEvents();
    
private:
    // File paths
    std::vector<std::string> kpsreco_files_;
    std::string truth_file_;
    
    // Configuration
    bool is_mc_;
    bool is_dlana_;
    
    // File managers
    std::unique_ptr<larlite::storage_manager> larlite_io_;
    std::unique_ptr<larcv::IOManager> larcv_io_;
    std::unique_ptr<TChain> kpsreco_;

    // Event tracking
    int current_entry_;
    int total_entries_;
    int kpsreco_run;
    int kpsreco_subrun;
    int kpsreco_event;
    
    // Helper methods
    bool setupLarliteIO();
    bool setupLarcvIO();
    bool setupRecoIO();
    bool parseTruthFileList();
    std::vector<std::string> truth_files_;
    
    // File validation
    bool validateFiles() const;
};

} // namespace gen2ntuple