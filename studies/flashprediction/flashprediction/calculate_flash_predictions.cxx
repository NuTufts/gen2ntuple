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

// ublarcvapp (for MC truth)
#include "ublarcvapp/MCTools/NeutrinoVertex.h"

// larlite (for MC truth data)
#include "DataFormat/mctruth.h"

// larutil (for SCE correction)
#include "LArUtil/SpaceChargeMicroBooNE.h"

// Flashmatch library to be able to run siren model and calculate sinkhorn divergences
// Torch headers
#include <torch/torch.h>

// Project headers
#include "flashmatch_dataprep/SirenTorchModel.h"
#include "flashmatch_dataprep/ModelInputInterface.h"
#include "flashmatch_dataprep/UnbalancedSinkhornDivergence.h"
#include "flashmatch_dataprep/UBFlashSinkDiv.h"

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
    std::cout << "  --siren-model-file <file>   Provide path to SIREN model file and activate SIREN flash prediction" << std::endl;
    std::cout << "  -h, --help                  Show this help message" << std::endl;
    
}

int main(int argc, char** argv) {
    
    // Parse command line arguments
    std::string dlmerged_file = "";
    std::string reco_file = "";
    std::string output_file = "";
    std::string siren_model_file = "";
    int num_entries = -1;
    int start_entry = 0;
    float adc_threshold = 10.0;
    bool verbose = false;
    bool tickbackward = false;
    bool is_mc = false;
    bool run_siren = false;
    float siren_pe_scale = 3.0;
    const int max_vertices = 5;
    
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
        else if (arg == "--siren-model-file") {
            run_siren = true;
            siren_model_file = argv[++i];
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

    if ( run_siren && siren_model_file.empty()) {
        std::cerr << "Error: SIREN model file provided but is empty." << std::endl;
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
    if ( run_siren ) {
        std::cout << "Siren model file: " << siren_model_file << std::endl;
    }
    else {
        std::cout << "Not running Siren Model" << std::endl;
    }
    
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

    // We ave an entry per event in this tree
    TTree* output_tree = new TTree("FlashPredictionTree", "Flash predictions for neutrino vertices");
    
    // --------------------------------------------------------------------
    // Variables to store in the TTree Branches

    // Basic event info
    int entry, run, subrun, event;
    int n_vertices;
    bool has_vertices, has_flash;
    
    // Observed flash info (same for whole event)
    float obs_total_pe, obs_time;
    std::vector<float> obs_pe_per_pmt;
    
    // Vectors for vertex-specific predictions (all particles)
    std::vector<float> reco_vertex_x_v; // reconstructed vertex x-position
    std::vector<float> reco_vertex_y_v; // reconstructed vertex y-position
    std::vector<float> reco_vertex_z_v; // reconstructed vertex z-position
    std::vector<int>   n_tracks_all_v;  // number of tracks in the vertex
    std::vector<int>   n_showers_all_v; // number of showers in the vertex
    std::vector<int>   n_primary_tracks_v;  // number of primary tracks in the vertex
    std::vector<int>   n_primary_showers_v; // number of primary showers in the vertex
    std::vector<float> total_charge_all_v;
    std::vector<float> total_photons_all_v;

    std::vector<float>              ubpred_total_pe_all_v;   // ub light model prediction, total pe
    std::vector<std::vector<float>> ubpred_pe_per_pmt_all_v; // ub light model pe per pmt prediction: [vertex][pmt]

    std::vector<float>              siren_total_pe_all_v; // ub light model prediction, total pe
    std::vector<std::vector<float>> siren_pe_per_pmt_all_v; // ub light model pe per pmt prediction: [vertex][pmt]
    
    // // Metrics vectors for UB light model
    std::vector<std::vector<float>> ub_sinkhorn_div_all_v;             // balanced sinkhorn divergance: [vertex][reg_param]
    std::vector<std::vector<float>> ub_unbalanced_sinkhorn_div_all_v;  // unbalanced sinkhorn divergence: [vertex][reg_param]
    std::vector<float>              ub_pe_diff_all_v;
    std::vector<float>              ub_pe_fracerr_all_v;

    // Metrics vectors for siren light model
    std::vector<std::vector<float>> siren_sinkhorn_div_all_v;            // balanced sinkhorn divergence: [vertex][reg_param]
    std::vector<std::vector<float>> siren_unbalanced_sinkhorn_div_all_v; // unbalanced sinkhorn divergence: [vertex][reg_param]
    std::vector<float>              siren_pe_diff_all_v;
    std::vector<float>              siren_pe_fracerr_all_v;
    
    // MC truth vectors (only used if is_mc is true)
    std::vector<float> vtx_dist_to_true_v;
    float true_vtx_x, true_vtx_y, true_vtx_z;
    bool has_mc_truth;
    
    // -----------------------------------------------------------------------------
    // Set up branches

    // Event-level branches
    output_tree->Branch("entry", &entry, "entry/I");
    output_tree->Branch("run", &run, "run/I");
    output_tree->Branch("subrun", &subrun, "subrun/I");
    output_tree->Branch("event", &event, "event/I");
    output_tree->Branch("n_vertices",   &n_vertices,   "n_vertices/I");
    output_tree->Branch("has_vertices", &has_vertices, "has_vertices/O");
    output_tree->Branch("has_flash",    &has_flash,    "has_flash/O");

    // Observed flash branches (event-level)
    output_tree->Branch("obs_total_pe",   &obs_total_pe,    "obs_total_pe/F");
    output_tree->Branch("obs_time",       &obs_time,        "obs_time/F");
    output_tree->Branch("obs_pe_per_pmt", &obs_pe_per_pmt);

    // Reconstructed vertex position branches (vectors)
    output_tree->Branch("reco_vertex_x", &reco_vertex_x_v);
    output_tree->Branch("reco_vertex_y", &reco_vertex_y_v);
    output_tree->Branch("reco_vertex_z", &reco_vertex_z_v);

    // Particle count branches for all particles (vectors)
    output_tree->Branch("n_tracks_all", &n_tracks_all_v);
    output_tree->Branch("n_showers_all", &n_showers_all_v);
    output_tree->Branch("n_primary_tracks", &n_primary_tracks_v);
    output_tree->Branch("n_primary_showers", &n_primary_showers_v);
    output_tree->Branch("total_charge_all", &total_charge_all_v);
    output_tree->Branch("total_photons_all", &total_photons_all_v);

    // UB light model prediction branches (vectors)
    output_tree->Branch("ubpred_total_pe_all", &ubpred_total_pe_all_v);
    output_tree->Branch("ubpred_pe_per_pmt_all", &ubpred_pe_per_pmt_all_v);

    // UB light model metrics branches (vectors)
    output_tree->Branch("ub_sinkhorn_div_all", &ub_sinkhorn_div_all_v);
    output_tree->Branch("ub_unbalanced_sinkhorn_div_all", &ub_unbalanced_sinkhorn_div_all_v);
    output_tree->Branch("ub_pe_diff_all", &ub_pe_diff_all_v);
    output_tree->Branch("ub_pe_fracerr_all", &ub_pe_fracerr_all_v);

    // SIREN model prediction branches (vectors)
    output_tree->Branch("siren_total_pe_all", &siren_total_pe_all_v);
    output_tree->Branch("siren_pe_per_pmt_all", &siren_pe_per_pmt_all_v);

    // SIREN model metrics branches (vectors)
    output_tree->Branch("siren_sinkhorn_div_all", &siren_sinkhorn_div_all_v);
    output_tree->Branch("siren_unbalanced_sinkhorn_div_all", &siren_unbalanced_sinkhorn_div_all_v);
    output_tree->Branch("siren_pe_diff_all", &siren_pe_diff_all_v);
    output_tree->Branch("siren_pe_fracerr_all", &siren_pe_fracerr_all_v);

    // MC truth branches (will fill with dummy values if no MC truth)
    output_tree->Branch("vtx_dist_to_true", &vtx_dist_to_true_v);
    output_tree->Branch("true_vtx_x", &true_vtx_x, "true_vtx_x/F");
    output_tree->Branch("true_vtx_y", &true_vtx_y, "true_vtx_y/F");
    output_tree->Branch("true_vtx_z", &true_vtx_z, "true_vtx_z/F");
    output_tree->Branch("has_mc_truth", &has_mc_truth, "has_mc_truth/O");


    // -----------------------------------------------------------------------------
    
    // --------------------------------------------------------------
    // UB light model and larflow::sinkhorn
    // Initialize flash predictor and Sinkhorn calculator
    larflow::reco::NuVertexFlashPrediction predictor;
    
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
    
    // Regularization parameters for Sinkhorn divergence (dampening)
    float sinkhorn_regularizations[3] = {0.1, 1.0, 10.0};

    // --------------------------------------------------------------
    // torch siren model and c++ implementation of geomloss sinkhorn

    // Create SirenTorchModel and load weights (only if run_siren is true)
    flashmatch::SirenTorchModel siren_model;
    flashmatch::ModelInputInterface input_interface;

    if (run_siren) {
        std::cout << "Loading Siren model from: " << siren_model_file << std::endl;
        if (verbose) {
            siren_model.set_verbosity(1);
        }
        try {
            siren_model.load_model_file(siren_model_file);
        }
        catch ( std::exception& e ) {
            std::cerr << "Could not load SIREN Model" << std::endl;
            std::cerr << e.what() << std::endl;
            return 0;
        }

        // Set normalization parameters for SIREN model input
        std::vector<float> planecharge_offset = {0.0, 0.0, 0.0};
        std::vector<float> planecharge_scale = {50000.0, 50000.0, 50000.0};
        input_interface.set_planecharge_normalization(planecharge_offset, planecharge_scale);
        input_interface.set_use_log_normalization(false);
    }

    // Create UBFlashSinkDiv for Sinkhorn divergence calculations
    flashmatch::UBFlashSinkDiv ubsinkdiv_algo;
    
    // --------------------------------------------------------------
    // Initialize MC truth tools (only if MC mode enabled)
    ublarcvapp::mctools::NeutrinoVertex* mc_nu_vertexer = nullptr;
    larutil::SpaceChargeMicroBooNE* sce = nullptr;  // goes from true position to space-charge modified position
    
    
    if (is_mc) {
        mc_nu_vertexer = new ublarcvapp::mctools::NeutrinoVertex();
        sce = new larutil::SpaceChargeMicroBooNE();
        
        if (verbose) {
            std::cout << "MC truth tools initialized" << std::endl;
        }
    }

    // Initialize Space-charge correction tool for reconstructed positions
    larutil::SpaceChargeMicroBooNE* reverse_sce = nullptr; // goes from observed position to space-charge corrected position
    reverse_sce = new larutil::SpaceChargeMicroBooNE( larutil::SpaceChargeMicroBooNE::kMCC9_Backward );
    
    // Process entries
    for (int ientry = start_entry; ientry < end_entry; ientry++) {
        
        if (verbose || ientry % 100 == 0) {
            std::cout << "Processing entry " << ientry << " / " << end_entry - 1 << std::endl;
        }
        
        // ------------------------------------------------
        // Clear all vectors for this event

        // Reconstructed vertex positions
        reco_vertex_x_v.clear();
        reco_vertex_y_v.clear();
        reco_vertex_z_v.clear();

        // Particle counts
        n_tracks_all_v.clear();
        n_showers_all_v.clear();
        n_primary_tracks_v.clear();
        n_primary_showers_v.clear();
        total_charge_all_v.clear();
        total_photons_all_v.clear();

        // UB light model predictions
        ubpred_total_pe_all_v.clear();
        ubpred_pe_per_pmt_all_v.clear();

        // UB light model metrics
        ub_sinkhorn_div_all_v.clear();
        ub_unbalanced_sinkhorn_div_all_v.clear();
        ub_pe_diff_all_v.clear();
        ub_pe_fracerr_all_v.clear();

        // SIREN model predictions
        siren_total_pe_all_v.clear();
        siren_pe_per_pmt_all_v.clear();

        // SIREN model metrics
        siren_sinkhorn_div_all_v.clear();
        siren_unbalanced_sinkhorn_div_all_v.clear();
        siren_pe_diff_all_v.clear();
        siren_pe_fracerr_all_v.clear();

        // MC truth
        vtx_dist_to_true_v.clear();

        // ------------------------------------------------
        
        
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
        
        const std::vector<larcv::Image2D>& adc_v = ev_img->as_vector();
        
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

                const larflow::reco::NuVertexCandidate& vertex_candidate = nuvetoed_v->at(vtx_idx);

                // Store the reconstructed vertex position
                reco_vertex_x_v.push_back(vertex_candidate.pos[0]);
                reco_vertex_y_v.push_back(vertex_candidate.pos[1]);
                reco_vertex_z_v.push_back(vertex_candidate.pos[2]);

                // Calculate distance to true vertex (only if MC mode enabled and MC truth available)
                if (is_mc && has_mc_truth) {
                    float dx = vertex_candidate.pos[0] - true_vtx_pos.X();
                    float dy = vertex_candidate.pos[1] - true_vtx_pos.Y();
                    float dz = vertex_candidate.pos[2] - true_vtx_pos.Z();
                    float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    vtx_dist_to_true_v.push_back(distance);

                    if (verbose) {
                        std::cout << "  Vertex[" << vtx_idx << "] distance to true: " << distance << " cm" << std::endl;
                    }
                } else if (is_mc) {
                    vtx_dist_to_true_v.push_back(-999.0);
                }

                // ========================================================================
                // UB Light Model Flash Prediction (using larflow predictor)
                // ========================================================================

                bool ubpred_success = false;
                std::vector<float> ubpred_pe_per_pmt(32, 0.0);
                float ubpred_total_pe = 0.0;

                try {
                    auto predicted_flash = predictor.predictFlash(
                        vertex_candidate,
                        adc_v,
                        adc_threshold,
                        true,   // use_trilinear
                        false   // primary_prongs_only = false (all particles)
                    );

                    ubpred_success = true;
                    ubpred_total_pe = predictor.getTotalPredictedPE();
                    n_tracks_all_v.push_back(predictor.getNumTracksProcessed());
                    n_showers_all_v.push_back(predictor.getNumShowersProcessed());
                    total_charge_all_v.push_back(predictor.getTotalChargeCollected());
                    total_photons_all_v.push_back(predictor.getTotalPhotonsEmitted());

                    // Count primary tracks and showers
                    int n_primary_tracks = 0;
                    int n_primary_showers = 0;
                    for (auto const& issecondary : vertex_candidate.track_isSecondary_v) {
                        if (issecondary==0) n_primary_tracks++;
                    }
                    for (auto const& issecondary : vertex_candidate.shower_isSecondary_v) {
                        if (issecondary==0) n_primary_showers++;
                    }
                    n_primary_tracks_v.push_back(n_primary_tracks);
                    n_primary_showers_v.push_back(n_primary_showers);

                    // Get per-PMT predictions
                    const auto& pe_per_pmt_map = predictor.getPredictedPE();
                    for (int pmt = 0; pmt < 32; pmt++) {
                        auto it = pe_per_pmt_map.find(pmt);
                        ubpred_pe_per_pmt[pmt] = (it != pe_per_pmt_map.end()) ? it->second : 0.0;
                    }

                } catch (const std::exception& e) {
                    if (verbose) {
                        std::cerr << "Warning: UB flash prediction failed for entry " << ientry
                                  << ", vertex " << vtx_idx << ": " << e.what() << std::endl;
                    }
                    n_tracks_all_v.push_back(0);
                    n_showers_all_v.push_back(0);
                    n_primary_tracks_v.push_back(0);
                    n_primary_showers_v.push_back(0);
                    total_charge_all_v.push_back(0.0);
                    total_photons_all_v.push_back(0.0);
                }

                ubpred_total_pe_all_v.push_back(ubpred_total_pe);
                ubpred_pe_per_pmt_all_v.push_back(ubpred_pe_per_pmt);

                // Calculate UB metrics
                float ub_pe_diff = ubpred_total_pe - obs_total_pe;
                float ub_pe_fracerr = (obs_total_pe > 0.0) ? (ubpred_total_pe - obs_total_pe) / obs_total_pe : -999.0;
                ub_pe_diff_all_v.push_back(ub_pe_diff);
                ub_pe_fracerr_all_v.push_back(ub_pe_fracerr);

                // Calculate UB Sinkhorn divergences (balanced and unbalanced)
                std::vector<float> ub_sinkhorn_balanced(1, -999.0);
                std::vector<float> ub_sinkhorn_unbalanced(1, -999.0);

                if (ubpred_success && has_flash) {
                    try {
                        ub_sinkhorn_balanced[0] = ubsinkdiv_algo.calc(ubpred_pe_per_pmt, obs_pe_per_pmt, true);
                    } catch (const std::exception& e) {
                        if (verbose) {
                            std::cerr << "Warning: UB balanced Sinkhorn failed: " << e.what() << std::endl;
                        }
                    }

                    try {
                        ub_sinkhorn_unbalanced[0] = ubsinkdiv_algo.calc(ubpred_pe_per_pmt, obs_pe_per_pmt, false);
                    } catch (const std::exception& e) {
                        if (verbose) {
                            std::cerr << "Warning: UB unbalanced Sinkhorn failed: " << e.what() << std::endl;
                        }
                    }
                }

                ub_sinkhorn_div_all_v.push_back(ub_sinkhorn_balanced);
                ub_unbalanced_sinkhorn_div_all_v.push_back(ub_sinkhorn_unbalanced);

                // ========================================================================
                // SIREN Model Flash Prediction
                // ========================================================================

                bool siren_success = false;
                std::vector<float> siren_pe_per_pmt(32, 0.0);
                float siren_total_pe = 0.0;

                if (run_siren) {
                    try {
                        // Get 3D points and charge from the vertex candidate's track and shower collections
                        std::vector<std::vector<float>> voxel_positions;
                        std::vector<std::vector<float>> voxel_pixelcoords;

                        // Extract hits from tracks and showers
                        std::vector< const std::vector<larlite::larflowcluster>* > phitclusters = {
                            &vertex_candidate.track_hitcluster_v,
                            &vertex_candidate.shower_v
                        };

                        for ( auto const& phitcluster : phitclusters ) {
                            for (auto const& trackcluster : *phitcluster ) {
                                for (size_t ipt = 0; ipt < trackcluster.size(); ipt++) {
                                    auto const& lfhit = trackcluster.at(ipt);

                                    // these are the 3D positions in the detector
                                    std::vector<float> pos = { lfhit[0], lfhit[1], lfhit[2] };

                                    // we have to space charge correct the positions
                                    bool applied = false;
                                    std::vector<double> pos_sce = reverse_sce->ApplySpaceChargeEffect( lfhit[0], lfhit[1], lfhit[2], applied );
                                    if ( !applied ) {
                                        // if a correction was not applied, this hit is not inside the TPC
                                        // we can throw it out.
                                        continue;
                                    }
                                    std::vector<float> fpos_sce(3,0);
                                    for (size_t v=0; v<3; v++)
                                        fpos_sce[v] = pos_sce[v];
                                    voxel_positions.push_back( fpos_sce );

                                    // these are the image coordinates from which they were projected
                                    std::vector<float> pixelcoords(4,0); // (tick, U,V,Y)
                                    pixelcoords[0] = lfhit.tick;
                                    pixelcoords[1] = lfhit.targetwire[0];
                                    pixelcoords[2] = lfhit.targetwire[1];
                                    pixelcoords[3] = lfhit.targetwire[2];

                                    voxel_pixelcoords.push_back(pixelcoords);

                                }
                            }
                        }

                        int num_voxels = voxel_positions.size();

                        if (num_voxels > 0) {
                            // Prepare input tensors for SIREN model
                            torch::Tensor voxel_features_t;
                            torch::Tensor voxel_charge_t; 

                            input_interface.prepare_input_tensor( 
                                voxel_positions, 
                                voxel_pixelcoords,
                                adc_v,
                                voxel_features_t,
                                voxel_charge_t
                            );

                            // Run SIREN model
                            std::vector<float> siren_output = siren_model.predict_pe(voxel_features_t, voxel_charge_t);
                            siren_pe_per_pmt.resize(32,0);
                            for (size_t ipmt=0; ipmt<32; ipmt++) {
                                siren_pe_per_pmt[ipmt] = siren_output[ipmt]*siren_pe_scale;
                            }

                            // Calculate total PE
                            for (int i = 0; i < 32; i++) {
                                siren_total_pe += siren_pe_per_pmt[i];
                            }

                            siren_success = true;
                        }// if num_voxels = 0

                    } catch (const std::exception& e) {
                        if (verbose) {
                            std::cerr << "Warning: SIREN prediction failed for entry " << ientry
                                      << ", vertex " << vtx_idx << ": " << e.what() << std::endl;
                        }
                    }

                }// if run_siren flag is True

                siren_total_pe_all_v.push_back(siren_total_pe);
                siren_pe_per_pmt_all_v.push_back(siren_pe_per_pmt);

                // Calculate SIREN metrics
                float siren_pe_diff = siren_total_pe - obs_total_pe;
                float siren_pe_fracerr = (obs_total_pe > 0.0) ? (siren_total_pe - obs_total_pe) / obs_total_pe : -999.0;
                siren_pe_diff_all_v.push_back(siren_pe_diff);
                siren_pe_fracerr_all_v.push_back(siren_pe_fracerr);

                // Calculate SIREN Sinkhorn divergences (balanced and unbalanced)
                std::vector<float> siren_sinkhorn_balanced(1, -999.0);
                std::vector<float> siren_sinkhorn_unbalanced(1, -999.0);

                if (siren_success && has_flash) {
                    try {
                        siren_sinkhorn_balanced[0] = ubsinkdiv_algo.calc(siren_pe_per_pmt, obs_pe_per_pmt, true);
                    } catch (const std::exception& e) {
                        if (verbose) {
                            std::cerr << "Warning: SIREN balanced Sinkhorn failed: " << e.what() << std::endl;
                        }
                    }

                    try {
                        siren_sinkhorn_unbalanced[0] = ubsinkdiv_algo.calc(siren_pe_per_pmt, obs_pe_per_pmt, false);
                    } catch (const std::exception& e) {
                        if (verbose) {
                            std::cerr << "Warning: SIREN unbalanced Sinkhorn failed: " << e.what() << std::endl;
                        }
                    }
                }

                siren_sinkhorn_div_all_v.push_back(siren_sinkhorn_balanced);
                siren_unbalanced_sinkhorn_div_all_v.push_back(siren_sinkhorn_unbalanced);

                // Verbose output
                if (verbose) {
                    std::cout << "  Vertex[" << vtx_idx << "]";
                    std::cout << " UB_PE=" << ubpred_total_pe;
                    if (run_siren) std::cout << " SIREN_PE=" << siren_total_pe;
                    std::cout << " OBS_PE=" << obs_total_pe;
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