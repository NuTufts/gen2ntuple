#include "FileManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>

// LArLite includes
#include "DataFormat/storage_manager.h"

// LArCV includes  
#include "larcv/core/DataFormat/IOManager.h"

namespace gen2ntuple {

FileManager::FileManager() 
    : is_mc_(false), 
    is_dlana_(false), 
    _nuvtx_v(nullptr),
    current_entry_(0), 
    total_entries_(0) {
}

FileManager::~FileManager() {
    closeFiles();
}

bool FileManager::setInputFiles(const std::vector<std::string>& kpsreco_files, 
                               const std::string& truth_file) {
    kpsreco_files_ = kpsreco_files;
    truth_file_ = truth_file;
    
    if (!validateFiles()) {
        return false;
    }
    
    if (!parseTruthFileList()) {
        return false;
    }
    
    return true;
}

bool FileManager::validateFiles() const {
    auto file_exists = [](const std::string& path) {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    };
    
    // Check KPSReco files
    for (const auto& file : kpsreco_files_) {
        if (!file_exists(file)) {
            std::cerr << "FileManager: KPSReco file does not exist: " << file << std::endl;
            return false;
        }
    }
    
    // Check truth file
    if (!file_exists(truth_file_)) {
        std::cerr << "FileManager: Truth file does not exist: " << truth_file_ << std::endl;
        return false;
    }
    
    return true;
}

bool FileManager::parseTruthFileList() {
    truth_files_.clear();
    
    // Check if truth_file_ is a list file (.txt) or a single ROOT file
    if (truth_file_.size() > 4 && truth_file_.substr(truth_file_.size() - 4) == ".txt") {
        // Parse list file
        std::ifstream file(truth_file_);
        if (!file.is_open()) {
            std::cerr << "FileManager: Cannot open truth file list: " << truth_file_ << std::endl;
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '#') {
                truth_files_.push_back(line);
            }
        }
        file.close();
    } else {
        // Single ROOT file
        truth_files_.push_back(truth_file_);
    }
    
    // Validate truth files
    auto file_exists = [](const std::string& path) {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    };
    
    for (const auto& file : truth_files_) {
        if (!file_exists(file)) {
            std::cerr << "FileManager: Truth file does not exist: " << file << std::endl;
            return false;
        }
    }
    
    return true;
}

bool FileManager::openFiles() {
    if (!setupLarliteIO()) {
        return false;
    }
    
    if (!setupLarcvIO()) {
        return false;
    }
    
    // Get total entries (minimum of both managers)
    int larlite_entries = larlite_io_->get_entries();
    int larcv_entries = larcv_io_->get_n_entries();
    
    total_entries_ = std::min(larlite_entries, larcv_entries);
    current_entry_ = 0;
    
    std::cout << "FileManager: Opened files successfully" << std::endl;
    std::cout << "  LArLite entries: " << larlite_entries << std::endl;
    std::cout << "  LArCV entries: " << larcv_entries << std::endl;
    std::cout << "  Processing entries: " << total_entries_ << std::endl;
    
    return true;
}

bool FileManager::setupLarliteIO() {
    larlite_io_ = std::make_unique<larlite::storage_manager>();
    
    // Set mode
    larlite_io_->set_io_mode(larlite::storage_manager::kREAD);
    
    // Add truth files (merged_dlreco)
    for (const auto& file : truth_files_) {
        larlite_io_->add_in_filename(file);
    }
    
    // Open
    if (!larlite_io_->open()) {
        std::cerr << "FileManager: Failed to open LArLite files" << std::endl;
        return false;
    }
    
    return true;
}

bool FileManager::setupLarcvIO() {
    auto tick_direction = larcv::IOManager::kTickBackward;

    larcv_io_ = std::make_unique<larcv::IOManager>(larcv::IOManager::kREAD,"",tick_direction);
    
    // Add truth files (merged_dlreco) - same files as LArLite, different trees
    for (const auto& file : truth_files_) {
        larcv_io_->add_in_file(file);
    }
    if (tick_direction==larcv::IOManager::kTickBackward) {
        larcv_io_->reverse_all_products();
    }
    
    // Initialize
    larcv_io_->initialize();
    
    return true;
}

bool FileManager::setupRecoIO() {

    kpsreco_ = std::make_unique<TChain>("KPSRecoManagerTree");
    
    // Add truth files (merged_dlreco)
    for (const auto& recofile : kpsreco_files_ ) {
        kpsreco_->Add( recofile.c_str() );
    }

    return true;
}

void FileManager::closeFiles() {
    if (larlite_io_) {
        larlite_io_->close();
        larlite_io_.reset();
    }
    
    if (larcv_io_) {
        larcv_io_->finalize();
        larcv_io_.reset();
    }

    if (kpsreco_) {
        kpsreco_.reset();
    }
}

bool FileManager::nextEvent() {
    if (isEndOfFile()) {
        return false;
    }
    
    // Read next event from both managers
    if (!larlite_io_->go_to(current_entry_)) {
        std::cerr << "FileManager: Failed to read LArLite entry " << current_entry_ << std::endl;
        return false;
    }
    
    larcv_io_->read_entry(current_entry_);

    kpsreco_->GetEntry(current_entry_);
    
    current_entry_++;
    
    // Verify event synchronization (with error handling)
    // Temporarily disabled to debug memory corruption
    /*
    try {
        if (!synchronizeEvents()) {
            std::cerr << "FileManager: Event synchronization failed at entry " 
                      << current_entry_ - 1 << std::endl;
            return false;
        }
    } catch (const std::exception& e) {
        std::cerr << "FileManager: Exception during event synchronization: " << e.what() << std::endl;
        return false;
    }
    */
    
    return true;
}

void FileManager::reset() {
    current_entry_ = 0;
}

FileManager::EventID FileManager::getCurrentEventID() const {
    EventID id;
    
    if (larlite_io_) {
        id.run = larlite_io_->run_id();
        id.subrun = larlite_io_->subrun_id();
        id.event = larlite_io_->event_id();
    } else {
        id.run = id.subrun = id.event = -1;
    }
    
    return id;
}

bool FileManager::synchronizeEvents() {

    if (!larlite_io_ || !larcv_io_ || kpsreco_ ) {
        return false;
    }
    
    try {
        // Get event IDs from both managers
        EventID larlite_id = getCurrentEventID();
        
        EventID larcv_id;
        larcv_id.run = larcv_io_->event_id().run();
        larcv_id.subrun = larcv_io_->event_id().subrun();  
        larcv_id.event = larcv_io_->event_id().event();

        EventID reco_id;
        reco_id.run    = kpsreco_run;
        reco_id.subrun = kpsreco_subrun;
        reco_id.event  = kpsreco_event;
        
        // Check if they match
        bool synchronized = (larlite_id == larcv_id) && ( larlite_id==reco_id );
        
        if (!synchronized) {
            std::cerr << "FileManager: Event ID mismatch!" << std::endl;
            std::cerr << "  LArLite: " << larlite_id.run << ":" << larlite_id.subrun 
                      << ":" << larlite_id.event << std::endl;
            std::cerr << "  LArCV: " << larcv_id.run << ":" << larcv_id.subrun 
                      << ":" << larcv_id.event << std::endl;
            std::cerr << "  Reco: " << reco_id.run << ":" << reco_id.subrun 
                      << ":" << reco_id.event << std::endl;
        }
        
        return synchronized;
    } catch (const std::exception& e) {
        std::cerr << "FileManager: Exception in synchronizeEvents: " << e.what() << std::endl;
        return false;
    }
}

} // namespace gen2ntuple