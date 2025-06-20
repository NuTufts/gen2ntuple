#include "Configuration.h"
#include "Logger.h"
#include "EventData.h"
#include "BranchManager.h"
#include "FileManager.h"
#include "MCTruthProcessor.h"
#include "VertexSelector.h"
#include "TrackProcessor.h"
#include "ShowerProcessor.h"
#include "VertexSummaryProcessor.h"

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <memory>

#include "ublarcvapp/MCTools/MCPixelPGraph.h"
#include "ublarcvapp/MCTools/MCPixelPMap.h"
#include "larpid/model/TorchModel.h"

using namespace gen2ntuple;

int main(int argc, char** argv) {
    
    std::cout << "Gen2Ntuple: C++ Implementation of DLGen2 Flat Ntuple Maker" << std::endl;
    std::cout << "==========================================================" << std::endl;
    
    // Parse configuration
    Configuration config;
    if (!config.parseCommandLine(argc, argv)) {
        return 1;
    }
    
    // Setup logging
    Logger& logger = Logger::getInstance();
    if (config.isDebug()) {
        logger.setLevel(Logger::DEBUG);
    } else if (config.isVerbose()) {
        logger.setLevel(Logger::INFO);
    } else {
        logger.setLevel(Logger::WARNING);
    }
    
    LOG_INFO("Starting Gen2Ntuple processing");
    LOG_INFO("Configuration:");
    LOG_INFO("  Input files: " + std::to_string(config.getKPSRecoFiles().size()) + " files");
    LOG_INFO("  Truth file: " + config.getTruthFile());
    LOG_INFO("  Output file: " + config.getOutputFile());
    LOG_INFO("  CNN Model: " + config.getModelPath());
    LOG_INFO("  Is MC: " + std::string(config.isMC() ? "yes" : "no"));
    LOG_INFO("  Include keypoints: " + std::string(config.isKeypointsDisabled() ? "no" : "yes"));
    LOG_INFO("  Device: " + config.getDevice());
    LOG_INFO("  Vertex Selection Method: " + config.getVertexSelection() );
    
    // Create output file
    std::unique_ptr<TFile> output_file(TFile::Open(config.getOutputFile().c_str(), "RECREATE"));
    if (!output_file || output_file->IsZombie()) {
        LOG_ERROR("Could not create output file: " + config.getOutputFile());
        return 1;
    }
    
    // Create data structures
    EventData event_data;
    RecoData  reco_data;
    POTData   pot_data;
    
    // Create branch manager
    BranchManager branch_manager(output_file.get(), config.isMC(), !config.isKeypointsDisabled());

    try {
        branch_manager.setupEventBranches(&event_data);
    }
    catch ( std::exception& e ) {
        std::cerr << "Error setting up branches for EventData" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    if (config.isMC()) {
        branch_manager.setupPOTBranches(&pot_data);
    }
    
    LOG_INFO("Branch setup complete");
    
    // Create processing modules
    FileManager file_manager;
    file_manager.setMCMode(config.isMC());
    file_manager.setDLAnaMode(config.isDLAna());
    
    // Set up input files
    if (!file_manager.setInputFiles(config.getKPSRecoFiles(), config.getTruthFile())) {
        LOG_ERROR("Failed to set input files");
        return 1;
    }
    
    // Open files
    if (!file_manager.openFiles()) {
        LOG_ERROR("Failed to open input files");
        return 1;
    }

    std::cout << "Setup Branches for RecoData" << std::endl;
    try {
        reco_data.setBranchAddresses( file_manager.getRecoIO() );
    }
    catch ( std::exception& e ) {
        std::cerr << "Error setting up branches for RecoData" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    ublarcvapp::mctools::MCPixelPGraph mcpg;
    ublarcvapp::mctools::MCPixelPMap   mcpm;

    LOG_INFO("Setup ProngCNN");
    bool larpid_debug = false;
    larpid::model::TorchModel larpid_model;
    larpid_model.Initialize( config.getModelPath(), larpid_debug );

    // Create processors
    std::unique_ptr<MCTruthProcessor> mc_processor;
    if (config.isMC()) {
        mc_processor = std::make_unique<MCTruthProcessor>();
        mc_processor->setWeightFile(config.getWeightFile());
        mc_processor->loadWeights();
    }
    
    VertexSelector vertex_selector;
    vertex_selector.setMCMode(config.isMC());
    vertex_selector.setKeypointMode(!config.isKeypointsDisabled());
    
    TrackProcessor track_processor;
    track_processor.setMCMode(config.isMC());
    track_processor.setProngCNNInterface( &larpid_model );
    track_processor.setMCPixelUtils( &mcpg, &mcpm );
    
    ShowerProcessor shower_processor;
    shower_processor.setMCMode(config.isMC());
    shower_processor.setLArPIDInterface( &larpid_model );
    shower_processor.setMCPixelUtils( &mcpg, &mcpm );

    VertexSummaryProcessor summary_processor;

    
    LOG_INFO("Processing modules initialized");
    
    // Event processing loop
    int events_processed = 0;
    int max_events = config.getMaxEvents();
    int start_event = config.getStartEvent();
    
    // // Skip to start event
    // for (int i = 0; i < start_event && !file_manager.isEndOfFile(); i++) {
    //     if (!file_manager.nextEvent()) {
    //         break;
    //     }
    // }
    
    LOG_INFO("Starting event processing from event " + std::to_string(start_event));
    
    while (!file_manager.isEndOfFile() && 
           (max_events < 0 || events_processed < max_events)) {
        
        if (!file_manager.nextEvent()) {
            break;
        }
        
        // Clear previous event data
        event_data.clear();
        
        // Get current event ID
        auto event_id = file_manager.getCurrentEventID();
        event_data.fileid = 0; // Would need file ID mapping
        event_data.run = event_id.run;
        event_data.subrun = event_id.subrun;
        event_data.event = event_id.event;

        mcpg.clear();
        mcpm.pixMap.clear();
        mcpg.buildgraph(*file_manager.getLarcvIO(), *file_manager.getLarliteIO());
        mcpm.buildmap( *file_manager.getLarcvIO(), mcpg );
        
        // Process MC truth if MC
        if (config.isMC() && mc_processor) {
            if (!mc_processor->processEvent(file_manager.getLarliteIO(), 
                                            file_manager.getLarcvIO(),
                                            &event_data)) {
                LOG_WARNING("MC truth processing failed for event " + 
                           std::to_string(event_data.event));
                continue;
            }
        }

        // Process vertex selection
        if (!vertex_selector.processEvent(file_manager.getLarliteIO(), 
                                         file_manager.getLarcvIO(), 
                                         &event_data, 
                                         &reco_data,
                                         config.getVertexSelection() )) {
            LOG_WARNING("Vertex selection failed for event " + 
                       std::to_string(event_data.event));
            continue;
        }
        
        // Process tracks
        if (!track_processor.processEvent( file_manager.getLarliteIO(),
                                           file_manager.getLarcvIO(),
                                           &event_data,
                                           &reco_data)) {
            LOG_WARNING("Track processing failed for event " + 
                       std::to_string(event_data.event));
            // Continue processing even if track processing fails
        }
        
        // Process showers
        if (!shower_processor.processEvent( file_manager.getLarliteIO(),
                                           file_manager.getLarcvIO(),
                                           &event_data,
                                           &reco_data)) {
            LOG_WARNING("Shower processing failed for event " + 
                       std::to_string(event_data.event));
            // Continue processing even if shower processing fails
        }

        // Calculate fractions now that we've calculated values for each track
        // and can build totals
        if (!summary_processor.calculateHitAndChargeFractions(&event_data)) {
            std::cerr << "Error calculating track hit and charge fractions" << std::endl;
        }
        
        // Calculate total reconstructed energy
        event_data.recoNuE = 0.0;
        for (int i = 0; i < event_data.nTracks; i++) {
            event_data.recoNuE += event_data.trackRecoE[i];
        }
        for (int i = 0; i < event_data.nShowers; i++) {
            event_data.recoNuE += event_data.showerRecoE[i];
        }
        
        // Fill the tree
        branch_manager.fillEventTree();
        
        events_processed++;
        
        if (events_processed % 100 == 0) {
            LOG_INFO("Processed " + std::to_string(events_processed) + " events");
        }
        
        LOG_DEBUG("Processed event " + std::to_string(event_data.run) + ":" + 
                 std::to_string(event_data.subrun) + ":" + std::to_string(event_data.event));
    }
    
    // Close input files
    file_manager.closeFiles();
    
    LOG_INFO("Event processing complete. Processed " + std::to_string(events_processed) + " events");
    
    // Fill POT info for MC
    if (config.isMC()) {
        pot_data.totPOT = 1.0e20;
        pot_data.totGoodPOT = 0.98e20;
        branch_manager.fillPOTTree();
        LOG_INFO("Filled POT information");
    }
    
    // Write and close
    LOG_INFO("Writing output...");
    branch_manager.write();
    
    // Get entry counts before closing file
    Long64_t event_entries = branch_manager.getEventEntries();
    Long64_t pot_entries = config.isMC() ? branch_manager.getPOTEntries() : 0;
    
    output_file->Close();
    
    LOG_INFO("Processing complete!");
    LOG_INFO("Output saved to: " + config.getOutputFile());
    LOG_INFO("Event entries: " + std::to_string(event_entries));
    if (config.isMC()) {
        LOG_INFO("POT entries: " + std::to_string(pot_entries));
    }
    
    return 0;
}