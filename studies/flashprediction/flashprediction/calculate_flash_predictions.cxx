/**
 * \file calculate_flash_predictions_vectorized.cxx
 *
 * \brief Executable to calculate flash predictions for neutrino vertex candidates
 *        with vectorized output to match ntuple structure
 *
 * This program:
 * 1. Loads dlmerged files (for ADC images and observed opflash)
 * 2. Loads reco analysis files (containing KPSRecoManagerTree and NuVertexCandidate objects)
 * 3. Calculates predicted flash for each neutrino vertex candidate
 * 4. Computes Sinkhorn divergence between predicted and observed flashes
 * 5. Saves results to a ROOT tree with vectors storing all vertices per event
 */

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

// larcv
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// larlite
#include "DataFormat/storage_manager.h"
#include "DataFormat/opflash.h"

// larflow
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/NuVertexFlashPrediction.h"
#include "larflow/Reco/SinkhornFlashDivergence.h"

// ublarcvapp (for MC truth)
#include "ublarcvapp/MCTools/NeutrinoVertex.h"

// larlite (for MC truth data)
#include "DataFormat/mctruth.h"

// larutil (for SCE correction)
#include "LArUtil/SpaceChargeMicroBooNE.h"

// ROOT (for TVector3)
#include "TVector3.h"

void printUsage() {
    std::cout << "Usage: calculate_flash_predictions [options]" << std::endl;
    std::cout << "\nRequired options:" << std::endl;
    std::cout << "  -d, --dlmerged     <file>    Input dlmerged ROOT file (larcv format)" << std::endl;
    std::cout << "  -r, --reco         <file>    Input reco file with KPSRecoManagerTree" << std::endl;
    std::cout << "  -o, --output       <file>    Output ROOT file for flash predictions" << std::endl;
    std::cout << "\nOptional:" << std::endl;
    std::cout << "  -n, --num-entries  <N>       Number of entries to process (default: all)" << std::endl;
    std::cout << "  -s, --start-entry  <N>       Starting entry (default: 0)" << std::endl;
    std::cout << "  -t, --threshold    <float>   ADC threshold (default: 10.0)" << std::endl;
    std::cout << "  -v, --verbose               Enable verbose output" << std::endl;
    std::cout << "  -tb, --tickbackward         Use tick backward direction" << std::endl;
    std::cout << "  -mc, --mc                   Enable MC mode (calculate distance to true vertex)" << std::endl;
    std::cout << "  -h, --help                  Show this help message" << std::endl;
}

int main(int argc, char** argv) {
    
    // Parse command line arguments
    std::string dlmerged_file = "";
    std::string reco_file = "";
    std::string output_file = "";
    int num_entries = -1;
    int start_entry = 0;
    float adc_threshold = 10.0;
    bool verbose = false;
    bool tickbackward = false;
    bool is_mc = false;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-d" || arg == "--dlmerged") {
            if (i + 1 < argc) dlmerged_file = argv[++i];
        }
        else if (arg == "-r" || arg == "--reco") {
            if (i + 1 < argc) reco_file = argv[++i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) output_file = argv[++i];
        }
        else if (arg == "-n" || arg == "--num-entries") {
            if (i + 1 < argc) num_entries = std::atoi(argv[++i]);
        }
        else if (arg == "-s" || arg == "--start-entry") {
            if (i + 1 < argc) start_entry = std::atoi(argv[++i]);
        }
        else if (arg == "-t" || arg == "--threshold") {
            if (i + 1 < argc) adc_threshold = std::atof(argv[++i]);
        }
        else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
        else if (arg == "-tb" || arg == "--tickbackward") {
            tickbackward = true;
        }
        else if (arg == "-mc" || arg == "--mc") {
            is_mc = true;
        }
        else if (arg == "-h" || arg == "--help") {
            printUsage();
            return 0;
        }
    }
    
    // Validate required arguments
    if (dlmerged_file.empty() || reco_file.empty() || output_file.empty()) {
        std::cerr << "Error: Missing required arguments!" << std::endl;
        printUsage();
        return 1;
    }
    
    std::cout << "Flash Prediction Calculator (Vectorized)" << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "DLMerged file: " << dlmerged_file << std::endl;
    std::cout << "Reco file: " << reco_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "ADC threshold: " << adc_threshold << std::endl;
    std::cout << "MC mode: " << (is_mc ? "enabled" : "disabled") << std::endl;
    
    // Open input files
    // 1. Open reco file with KPSRecoManagerTree
    TFile* reco_tfile = TFile::Open(reco_file.c_str(), "READ");
    if (!reco_tfile || reco_tfile->IsZombie()) {
        std::cerr << "Error: Cannot open reco file: " << reco_file << std::endl;
        return 1;
    }
    
    TTree* kps_tree = (TTree*)reco_tfile->Get("KPSRecoManagerTree");
    if (!kps_tree) {
        std::cerr << "Error: Cannot find KPSRecoManagerTree in " << reco_file << std::endl;
        return 1;
    }
    
    // Set up branch for vertex candidates
    std::vector<larflow::reco::NuVertexCandidate>* nuvetoed_v = nullptr;
    kps_tree->SetBranchAddress("nuvetoed_v", &nuvetoed_v);
    
    // 2. Set up larcv IOManager for ADC images
    auto tick_direction = larcv::IOManager::kTickForward;
    if (tickbackward)
        tick_direction = larcv::IOManager::kTickBackward;
    larcv::IOManager ioman(larcv::IOManager::kREAD, "ioman", tick_direction);
    ioman.add_in_file(dlmerged_file);
    if (tickbackward)
        ioman.reverse_all_products();
    ioman.initialize();
    
    // 3. Set up larlite storage manager for opflash
    larlite::storage_manager ioll(larlite::storage_manager::kREAD);
    ioll.add_in_filename(dlmerged_file);
    if (!ioll.open()) {
        std::cerr << "Error: Cannot open dlmerged file for larlite: " << dlmerged_file << std::endl;
        return 1;
    }
    
    // Determine number of entries to process
    int total_entries = kps_tree->GetEntries();
    int end_entry = (num_entries < 0) ? total_entries : std::min(start_entry + num_entries, total_entries);
    
    std::cout << "Processing entries " << start_entry << " to " << end_entry - 1 
              << " (total: " << end_entry - start_entry << ")" << std::endl;
    if ( is_mc ) {
        std::cout << "The file is indicated to be a simulated (aka a MC) file." << std::endl;
    }
    
    // Create output file and tree
    TFile* output_tfile = TFile::Open(output_file.c_str(), "RECREATE");
    TTree* output_tree = new TTree("FlashPredictionTree", "Flash predictions for neutrino vertices");
    
    // Define vector branches to store all vertices per event
    // Basic event info
    int entry, run, subrun, event;
    int n_vertices;
    bool has_vertices, has_flash;
    
    // Observed flash info (same for whole event)
    float obs_total_pe, obs_time;
    std::vector<float> obs_pe_per_pmt;
    
    // Vectors for vertex-specific predictions (all particles)
    std::vector<float> pred_total_pe_all_v;
    std::vector<std::vector<float>> pred_pe_per_pmt_all_v; // 2D vector: [vertex][pmt]
    std::vector<int> n_tracks_all_v;
    std::vector<int> n_showers_all_v;
    std::vector<float> total_charge_all_v;
    std::vector<float> total_photons_all_v;
    
    // Vectors for vertex-specific predictions (primary particles only)
    std::vector<float> pred_total_pe_primary_v;
    std::vector<std::vector<float>> pred_pe_per_pmt_primary_v; // 2D vector: [vertex][pmt]
    std::vector<int> n_tracks_primary_v;
    std::vector<int> n_showers_primary_v;
    std::vector<float> total_charge_primary_v;
    std::vector<float> total_photons_primary_v;
    
    // Metrics vectors
    std::vector<std::vector<float>> sinkhorn_div_all_v; // 2D vector: [vertex][reg_param]
    std::vector<std::vector<float>> sinkhorn_div_primary_v; // 2D vector: [vertex][reg_param]
    std::vector<float> pe_diff_all_v;
    std::vector<float> pe_diff_primary_v;
    std::vector<float> pe_ratio_all_v;
    std::vector<float> pe_ratio_primary_v;
    
    // Status vectors
    std::vector<bool> prediction_success_all_v;
    std::vector<bool> prediction_success_primary_v;
    
    // MC truth vectors (only used if is_mc is true)
    std::vector<float> vtx_dist_to_true_v;
    float true_vtx_x, true_vtx_y, true_vtx_z;
    bool has_mc_truth;
    
    // Set up branches
    // Event-level branches
    output_tree->Branch("entry", &entry, "entry/I");
    output_tree->Branch("run", &run, "run/I");
    output_tree->Branch("subrun", &subrun, "subrun/I");
    output_tree->Branch("event", &event, "event/I");
    output_tree->Branch("n_vertices", &n_vertices, "n_vertices/I");
    output_tree->Branch("has_vertices", &has_vertices, "has_vertices/O");
    output_tree->Branch("has_flash", &has_flash, "has_flash/O");
    
    // Observed flash branches (event-level)
    output_tree->Branch("obs_total_pe", &obs_total_pe, "obs_total_pe/F");
    output_tree->Branch("obs_time", &obs_time, "obs_time/F");
    output_tree->Branch("obs_pe_per_pmt", &obs_pe_per_pmt);
    
    // Prediction branches (vectors)
    output_tree->Branch("pred_total_pe_all", &pred_total_pe_all_v);
    output_tree->Branch("pred_pe_per_pmt_all", &pred_pe_per_pmt_all_v);
    output_tree->Branch("n_tracks_all", &n_tracks_all_v);
    output_tree->Branch("n_showers_all", &n_showers_all_v);
    output_tree->Branch("total_charge_all", &total_charge_all_v);
    output_tree->Branch("total_photons_all", &total_photons_all_v);
    
    output_tree->Branch("pred_total_pe_primary", &pred_total_pe_primary_v);
    output_tree->Branch("pred_pe_per_pmt_primary", &pred_pe_per_pmt_primary_v);
    output_tree->Branch("n_tracks_primary", &n_tracks_primary_v);
    output_tree->Branch("n_showers_primary", &n_showers_primary_v);
    output_tree->Branch("total_charge_primary", &total_charge_primary_v);
    output_tree->Branch("total_photons_primary", &total_photons_primary_v);
    
    // Metrics branches (vectors)
    output_tree->Branch("sinkhorn_div_all", &sinkhorn_div_all_v);
    output_tree->Branch("sinkhorn_div_primary", &sinkhorn_div_primary_v);
    output_tree->Branch("pe_diff_all", &pe_diff_all_v);
    output_tree->Branch("pe_diff_primary", &pe_diff_primary_v);
    output_tree->Branch("pe_ratio_all", &pe_ratio_all_v);
    output_tree->Branch("pe_ratio_primary", &pe_ratio_primary_v);
    
    // Status branches (vectors)
    output_tree->Branch("prediction_success_all", &prediction_success_all_v);
    output_tree->Branch("prediction_success_primary", &prediction_success_primary_v);
    
    // MC truth branches (only if MC mode enabled)
    if (is_mc) {
        output_tree->Branch("vtx_dist_to_true", &vtx_dist_to_true_v);
        output_tree->Branch("true_vtx_x", &true_vtx_x, "true_vtx_x/F");
        output_tree->Branch("true_vtx_y", &true_vtx_y, "true_vtx_y/F");
        output_tree->Branch("true_vtx_z", &true_vtx_z, "true_vtx_z/F");
        output_tree->Branch("has_mc_truth", &has_mc_truth, "has_mc_truth/O");
    }
    
    // Initialize flash predictor and Sinkhorn calculator
    larflow::reco::NuVertexFlashPrediction predictor;
    larflow::reco::SinkhornFlashDivergence sinkhorn_calc;
    //if (verbose)
    //   sinkhorn_calc.set_verbosity(larcv::msg::kINFO);
    
    // Configure flash predictor with standard parameters
    predictor.setChargeToPhotonParams(
        200.0,    // adc_per_electron
        23.6e-3,  // mev_per_electron (MeV)
        24000.0,  // photons_per_mev
        0.7       // recombination_factor
    );
    
    predictor.setTrackConversionParams(
        3,      // dcol
        3,      // drow
        0.3,    // minstepsize (cm)
        0.5     // maxstepsize (cm)
    );
    
    predictor.setShowerConversionParams(
        3,      // dcol
        3       // drow
    );
    
    // Regularization parameters for Sinkhorn divergence
    float sinkhorn_regularizations[3] = {0.1, 1.0, 10.0};
    
    // Initialize MC truth tools (only if MC mode enabled)
    ublarcvapp::mctools::NeutrinoVertex* mc_nu_vertexer = nullptr;
    larutil::SpaceChargeMicroBooNE* sce = nullptr;
    
    if (is_mc) {
        mc_nu_vertexer = new ublarcvapp::mctools::NeutrinoVertex();
        sce = new larutil::SpaceChargeMicroBooNE();
        
        if (verbose) {
            std::cout << "MC truth tools initialized" << std::endl;
        }
    }
    
    // Process entries
    for (int ientry = start_entry; ientry < end_entry; ientry++) {
        
        if (verbose || ientry % 100 == 0) {
            std::cout << "Processing entry " << ientry << " / " << end_entry - 1 << std::endl;
        }
        
        // Clear all vectors for this event
        pred_total_pe_all_v.clear();
        pred_pe_per_pmt_all_v.clear();
        n_tracks_all_v.clear();
        n_showers_all_v.clear();
        total_charge_all_v.clear();
        total_photons_all_v.clear();
        
        pred_total_pe_primary_v.clear();
        pred_pe_per_pmt_primary_v.clear();
        n_tracks_primary_v.clear();
        n_showers_primary_v.clear();
        total_charge_primary_v.clear();
        total_photons_primary_v.clear();
        
        sinkhorn_div_all_v.clear();
        sinkhorn_div_primary_v.clear();
        pe_diff_all_v.clear();
        pe_diff_primary_v.clear();
        pe_ratio_all_v.clear();
        pe_ratio_primary_v.clear();
        
        prediction_success_all_v.clear();
        prediction_success_primary_v.clear();
        
        // Clear MC truth vectors (only if MC mode enabled)
        if (is_mc) {
            vtx_dist_to_true_v.clear();
        }
        
        obs_pe_per_pmt.clear();
        obs_pe_per_pmt.resize(32, 0.0);
        
        // Set event info
        entry = ientry;
        
        // Read data from files
        kps_tree->GetEntry(ientry);
        ioman.read_entry(ientry);
        ioll.go_to(ientry);
        
        // Set run/subrun/event from larlite
        run = ioll.run_id();
        subrun = ioll.subrun_id();
        event = ioll.event_id();
        
        // Get MC truth information (only if MC mode enabled)
        TVector3 true_vtx_pos(0, 0, 0);
        has_mc_truth = false;
        
        if (is_mc) {
            try {
                // Get MC truth data
                auto ev_mctruth = (larlite::event_mctruth*)(ioll.get_data(larlite::data::kMCTruth, "generator"));
                std::cout << "ev_mctruth->size()=" << ev_mctruth->size() << std::endl;
                
                if (ev_mctruth && ev_mctruth->size() > 0) {
                    // Get true neutrino vertex position with SCE correction
                    std::vector<float> mc_nu_vertex(3,0);
                    mc_nu_vertex = mc_nu_vertexer->getPos3DwSCE(ioll, sce); // returns (x,y,z,tick)
                    
                    if (mc_nu_vertex.size() >= 3) {
                        true_vtx_pos.SetXYZ(mc_nu_vertex[0], mc_nu_vertex[1], mc_nu_vertex[2]);
                        true_vtx_x = true_vtx_pos.X();
                        true_vtx_y = true_vtx_pos.Y();
                        true_vtx_z = true_vtx_pos.Z();
                        has_mc_truth = true;
                        
                        if (verbose) {
                            std::cout << "True vertex at (" << true_vtx_x << ", " << true_vtx_y << ", " << true_vtx_z << ")" << std::endl;
                        }
                    }
                    else {
                        std::cerr << "Error: getPos3DwSCE returns a bad mc nu vertex position." << std::endl;
                        std::cerr << "mc_nu_vertex.size()==" << mc_nu_vertex.size() << std::endl;
                    }
                }
            } catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "Warning: MC truth processing failed for entry " << ientry << ": " << e.what() << std::endl;
                }
                has_mc_truth = false;
                true_vtx_x = -999.0;
                true_vtx_y = -999.0;
                true_vtx_z = -999.0;
            }
        }
        
        // Get ADC images
        auto ev_img = (larcv::EventImage2D*)(ioman.get_data(larcv::kProductImage2D, "wire"));
        if (!ev_img || ev_img->Image2DArray().size() < 3) {
            std::cerr << "Warning: Cannot get ADC images for entry " << ientry << std::endl;
            // Fill with defaults and continue
            n_vertices = 0;
            has_vertices = false;
            has_flash = false;
            obs_total_pe = 0.0;
            obs_time = -999.0;
            output_tree->Fill();
            continue;
        }
        
        std::vector<larcv::Image2D> adc_v;
        for (size_t p = 0; p < 3; p++) {
            adc_v.push_back(ev_img->Image2DArray()[p]);
        }
        
        // Get observed opflash
        auto ev_opflash = (larlite::event_opflash*)(ioll.get_data(larlite::data::kOpFlash, "simpleFlashBeam"));
        
        has_flash = (ev_opflash && ev_opflash->size() > 0);
        
        if (has_flash) {
            // Use the first flash (highest PE)
            const auto& flash = ev_opflash->at(0);
            obs_total_pe = flash.TotalPE();
            obs_time = flash.Time();
            
            for (int pmt = 0; pmt < 32; pmt++) {
                obs_pe_per_pmt[pmt] = flash.PE(pmt);
            }
        } else {
            // Create flat dummy opflash
            obs_total_pe = 0.0;
            obs_time = -1.0;
            for (int pmt = 0; pmt < 32; pmt++) {
                obs_pe_per_pmt[pmt] = 1.0/32.0;
            }
        }
        
        // Process vertex candidates
        has_vertices = (nuvetoed_v && nuvetoed_v->size() > 0);
        n_vertices = has_vertices ? nuvetoed_v->size() : 0;
        
        if (verbose) {
            if (has_vertices) {
                std::cout << "Entry " << ientry << ": " << n_vertices << " neutrino candidates" << std::endl;
            } else {
                std::cout << "Entry " << ientry << ": No neutrino candidates" << std::endl;
            }
        }
        
        // Initialize MC truth event-level variables when not in MC mode
        if (!is_mc) {
            has_mc_truth = false;
            true_vtx_x = -999.0;
            true_vtx_y = -999.0;
            true_vtx_z = -999.0;
        }
        
        if (has_vertices) {
            // Process each vertex candidate
            for (size_t vtx_idx = 0; vtx_idx < nuvetoed_v->size(); vtx_idx++) {
                
                const auto& vertex_candidate = nuvetoed_v->at(vtx_idx);
                
                // Initialize default values for this vertex
                pred_total_pe_all_v.push_back(-1.0);
                pred_total_pe_primary_v.push_back(-1.0);
                n_tracks_all_v.push_back(0);
                n_showers_all_v.push_back(0);
                n_tracks_primary_v.push_back(0);
                n_showers_primary_v.push_back(0);
                total_charge_all_v.push_back(0.0);
                total_photons_all_v.push_back(0.0);
                total_charge_primary_v.push_back(0.0);
                total_photons_primary_v.push_back(0.0);
                
                prediction_success_all_v.push_back(false);
                prediction_success_primary_v.push_back(false);
                
                // Initialize PMT vectors for this vertex
                std::vector<float> pmt_pe_all(32, 0.0);
                std::vector<float> pmt_pe_primary(32, 0.0);
                pred_pe_per_pmt_all_v.push_back(pmt_pe_all);
                pred_pe_per_pmt_primary_v.push_back(pmt_pe_primary);
                
                // Initialize metric vectors for this vertex
                std::vector<float> sinkhorn_all(3, -999.0);
                std::vector<float> sinkhorn_primary(3, -999.0);
                sinkhorn_div_all_v.push_back(sinkhorn_all);
                sinkhorn_div_primary_v.push_back(sinkhorn_primary);
                
                pe_diff_all_v.push_back(-999.0);
                pe_diff_primary_v.push_back(-999.0);
                pe_ratio_all_v.push_back(-999.0);
                pe_ratio_primary_v.push_back(-999.0);
                
                // Calculate distance to true vertex (only if MC mode enabled and MC truth available)
                if (is_mc && has_mc_truth) {
                    // Calculate 3D Euclidean distance
                    float dx = vertex_candidate.pos[0] - true_vtx_pos.X();
                    float dy = vertex_candidate.pos[1] - true_vtx_pos.Y();
                    float dz = vertex_candidate.pos[2] - true_vtx_pos.Z();
                    float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    vtx_dist_to_true_v.push_back(distance);
                    
                    if (verbose) {
                        std::cout << "  Vertex[" << vtx_idx << "] distance to true: " << distance << " cm" << std::endl;
                    }
                } else if (is_mc) {
                    // MC mode but no truth available
                    vtx_dist_to_true_v.push_back(-999.0);
                }
                
                // Prediction with all particles
                try {
                    auto predicted_flash_all = predictor.predictFlash(
                        vertex_candidate,
                        adc_v,
                        adc_threshold,
                        true,   // use_trilinear
                        false   // primary_prongs_only = false (all particles)
                    );
                    
                    prediction_success_all_v[vtx_idx] = true;
                    pred_total_pe_all_v[vtx_idx] = predictor.getTotalPredictedPE();
                    n_tracks_all_v[vtx_idx] = predictor.getNumTracksProcessed();
                    n_showers_all_v[vtx_idx] = predictor.getNumShowersProcessed();
                    total_charge_all_v[vtx_idx] = predictor.getTotalChargeCollected();
                    total_photons_all_v[vtx_idx] = predictor.getTotalPhotonsEmitted();
                    
                    // Get per-PMT predictions
                    const auto& pe_per_pmt_all = predictor.getPredictedPE();
                    for (int pmt = 0; pmt < 32; pmt++) {
                        auto it = pe_per_pmt_all.find(pmt);
                        pred_pe_per_pmt_all_v[vtx_idx][pmt] = (it != pe_per_pmt_all.end()) ? it->second : 0.0;
                    }
                    
                } catch (const std::exception& e) {
                    if (verbose) {
                        std::cerr << "Warning: Flash prediction (all) failed for entry " << ientry 
                                  << ", vertex " << vtx_idx << ": " << e.what() << std::endl;
                    }
                }
                
                // Prediction with primary particles only
                try {
                    auto predicted_flash_primary = predictor.predictFlash(
                        vertex_candidate,
                        adc_v,
                        adc_threshold,
                        true,   // use_trilinear
                        true    // primary_prongs_only = true
                    );
                    
                    prediction_success_primary_v[vtx_idx] = true;
                    pred_total_pe_primary_v[vtx_idx] = predictor.getTotalPredictedPE();
                    n_tracks_primary_v[vtx_idx] = predictor.getNumTracksProcessed();
                    n_showers_primary_v[vtx_idx] = predictor.getNumShowersProcessed();
                    total_charge_primary_v[vtx_idx] = predictor.getTotalChargeCollected();
                    total_photons_primary_v[vtx_idx] = predictor.getTotalPhotonsEmitted();
                    
                    // Get per-PMT predictions
                    const auto& pe_per_pmt_primary = predictor.getPredictedPE();
                    for (int pmt = 0; pmt < 32; pmt++) {
                        auto it = pe_per_pmt_primary.find(pmt);
                        pred_pe_per_pmt_primary_v[vtx_idx][pmt] = (it != pe_per_pmt_primary.end()) ? it->second : 0.0;
                    }
                    
                } catch (const std::exception& e) {
                    if (verbose) {
                        std::cerr << "Warning: Flash prediction (primary) failed for entry " << ientry 
                                  << ", vertex " << vtx_idx << ": " << e.what() << std::endl;
                    }
                }
                
                // Calculate metrics
                if (prediction_success_all_v[vtx_idx]) {
                    pe_diff_all_v[vtx_idx] = pred_total_pe_all_v[vtx_idx] - obs_total_pe;
                    if (obs_total_pe > 0.0) {
                        pe_ratio_all_v[vtx_idx] = pred_total_pe_all_v[vtx_idx] / obs_total_pe;
                    } else {
                        pe_ratio_all_v[vtx_idx] = 0.0;
                    }
                    
                    // Calculate Sinkhorn divergences for all particles
                    for (int i = 0; i < 3; i++) {
                        try {
                            sinkhorn_div_all_v[vtx_idx][i] = sinkhorn_calc.calculateDivergence(
                                pred_pe_per_pmt_all_v[vtx_idx],
                                obs_pe_per_pmt,
                                sinkhorn_regularizations[i],
                                100,    // max_iterations
                                1e-6    // tolerance
                            );
                        } catch (const std::exception& e) {
                            if (verbose) {
                                std::cerr << "Warning: Sinkhorn calculation failed (all, reg=" 
                                          << sinkhorn_regularizations[i] << "): " << e.what() << std::endl;
                            }
                            sinkhorn_div_all_v[vtx_idx][i] = -1.0; // Invalid value
                        }
                    }
                }
                
                if (prediction_success_primary_v[vtx_idx]) {
                    pe_diff_primary_v[vtx_idx] = pred_total_pe_primary_v[vtx_idx] - obs_total_pe;
                    if (obs_total_pe > 0.0) {
                        pe_ratio_primary_v[vtx_idx] = pred_total_pe_primary_v[vtx_idx] / obs_total_pe;
                    } else {
                        pe_ratio_primary_v[vtx_idx] = 0.0;
                    }
                    
                    // Calculate Sinkhorn divergences for primary particles
                    for (int i = 0; i < 3; i++) {
                        try {
                            sinkhorn_div_primary_v[vtx_idx][i] = sinkhorn_calc.calculateDivergence(
                                pred_pe_per_pmt_primary_v[vtx_idx],
                                obs_pe_per_pmt,
                                sinkhorn_regularizations[i],
                                100,    // max_iterations
                                1e-6    // tolerance
                            );
                        } catch (const std::exception& e) {
                            if (verbose) {
                                std::cerr << "Warning: Sinkhorn calculation failed (primary, reg=" 
                                          << sinkhorn_regularizations[i] << "): " << e.what() << std::endl;
                            }
                            sinkhorn_div_primary_v[vtx_idx][i] = -1.0; // Invalid value
                        }
                    }
                }
                
                if (verbose) {
                    std::cout << "  Vertex[" << vtx_idx << "] pred_PE(all)=" << pred_total_pe_all_v[vtx_idx] 
                              << " pred_PE(primary)=" << pred_total_pe_primary_v[vtx_idx];
                    if (is_mc && has_mc_truth) {
                        std::cout << " dist_to_true=" << vtx_dist_to_true_v[vtx_idx] << " cm";
                    }
                    std::cout << std::endl;
                }
                
            } // end loop over vertices
        } // end if has_vertices
        
        // Fill tree once per event
        output_tree->Fill();
        
    } // end loop over entries
    
    // Write output and cleanup
    output_tfile->cd();
    output_tree->Write();
    
    std::cout << "\nProcessing complete!" << std::endl;
    std::cout << "Output entries written: " << output_tree->GetEntries() << std::endl;
    
    output_tfile->Close();
    reco_tfile->Close();
    ioman.finalize();
    ioll.close();
    
    // Cleanup MC truth tools
    if (is_mc) {
        delete mc_nu_vertexer;
        delete sce;
    }
    
    std::cout << "\nOutput saved to: " << output_file << std::endl;
    
    return 0;
}