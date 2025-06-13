/**
 * \file calculate_flash_predictions.cxx
 *
 * \brief Executable to calculate flash predictions for neutrino vertex candidates
 *
 * This program:
 * 1. Loads dlmerged files (for ADC images and observed opflash)
 * 2. Loads reco analysis files (containing KPSRecoManagerTree and NuVertexCandidate objects)
 * 3. Calculates predicted flash for each neutrino vertex candidate
 * 4. Computes Sinkhorn divergence between predicted and observed flashes
 * 5. Saves results to a ROOT tree for friend tree analysis
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

struct FlashPredictionResult {
    int entry;
    int vertex_idx;
    
    // Predictions using all particles
    float pred_total_pe_all;
    float pred_pe_per_pmt_all[32];
    int n_tracks_all;
    int n_showers_all;
    float total_charge_all;
    float total_photons_all;
    
    // Predictions using primary particles only
    float pred_total_pe_primary;
    float pred_pe_per_pmt_primary[32];
    int n_tracks_primary;
    int n_showers_primary;
    float total_charge_primary;
    float total_photons_primary;
    
    // Observed flash info
    float obs_total_pe;
    float obs_pe_per_pmt[32];
    float obs_time;
    
    // Metrics
    float sinkhorn_div_all[3];      // For regularization params 0.1, 1.0, 10.0
    float sinkhorn_div_primary[3];
    float pe_diff_all;              // pred_total_pe_all - obs_total_pe
    float pe_diff_primary;
    float pe_ratio_all;             // pred_total_pe_all / obs_total_pe
    float pe_ratio_primary;
    
    // Status flags
    bool has_vertex;
    bool has_flash;
    bool prediction_success_all;
    bool prediction_success_primary;
};

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
    
    std::cout << "Flash Prediction Calculator" << std::endl;
    std::cout << "===========================" << std::endl;
    std::cout << "DLMerged file: " << dlmerged_file << std::endl;
    std::cout << "Reco file: " << reco_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "ADC threshold: " << adc_threshold << std::endl;
    
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
    larcv::IOManager ioman(larcv::IOManager::kREAD);
    ioman.add_in_file(dlmerged_file);
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
    
    // Create output file and tree
    TFile* output_tfile = TFile::Open(output_file.c_str(), "RECREATE");
    TTree* output_tree = new TTree("FlashPredictionTree", "Flash predictions for neutrino vertices");
    
    // Create result structure and set up branches
    FlashPredictionResult result;
    
    // Basic info branches
    output_tree->Branch("entry", &result.entry, "entry/I");
    output_tree->Branch("vertex_idx", &result.vertex_idx, "vertex_idx/I");
    
    // All particles prediction branches
    output_tree->Branch("pred_total_pe_all", &result.pred_total_pe_all, "pred_total_pe_all/F");
    output_tree->Branch("pred_pe_per_pmt_all", result.pred_pe_per_pmt_all, "pred_pe_per_pmt_all[32]/F");
    output_tree->Branch("n_tracks_all", &result.n_tracks_all, "n_tracks_all/I");
    output_tree->Branch("n_showers_all", &result.n_showers_all, "n_showers_all/I");
    output_tree->Branch("total_charge_all", &result.total_charge_all, "total_charge_all/F");
    output_tree->Branch("total_photons_all", &result.total_photons_all, "total_photons_all/F");
    
    // Primary particles prediction branches
    output_tree->Branch("pred_total_pe_primary", &result.pred_total_pe_primary, "pred_total_pe_primary/F");
    output_tree->Branch("pred_pe_per_pmt_primary", result.pred_pe_per_pmt_primary, "pred_pe_per_pmt_primary[32]/F");
    output_tree->Branch("n_tracks_primary", &result.n_tracks_primary, "n_tracks_primary/I");
    output_tree->Branch("n_showers_primary", &result.n_showers_primary, "n_showers_primary/I");
    output_tree->Branch("total_charge_primary", &result.total_charge_primary, "total_charge_primary/F");
    output_tree->Branch("total_photons_primary", &result.total_photons_primary, "total_photons_primary/F");
    
    // Observed flash branches
    output_tree->Branch("obs_total_pe", &result.obs_total_pe, "obs_total_pe/F");
    output_tree->Branch("obs_pe_per_pmt", result.obs_pe_per_pmt, "obs_pe_per_pmt[32]/F");
    output_tree->Branch("obs_time", &result.obs_time, "obs_time/F");
    
    // Metric branches
    output_tree->Branch("sinkhorn_div_all", result.sinkhorn_div_all, "sinkhorn_div_all[3]/F");
    output_tree->Branch("sinkhorn_div_primary", result.sinkhorn_div_primary, "sinkhorn_div_primary[3]/F");
    output_tree->Branch("pe_diff_all", &result.pe_diff_all, "pe_diff_all/F");
    output_tree->Branch("pe_diff_primary", &result.pe_diff_primary, "pe_diff_primary/F");
    output_tree->Branch("pe_ratio_all", &result.pe_ratio_all, "pe_ratio_all/F");
    output_tree->Branch("pe_ratio_primary", &result.pe_ratio_primary, "pe_ratio_primary/F");
    
    // Status branches
    output_tree->Branch("has_vertex", &result.has_vertex, "has_vertex/O");
    output_tree->Branch("has_flash", &result.has_flash, "has_flash/O");
    output_tree->Branch("prediction_success_all", &result.prediction_success_all, "prediction_success_all/O");
    output_tree->Branch("prediction_success_primary", &result.prediction_success_primary, "prediction_success_primary/O");
    
    // Initialize flash predictor and Sinkhorn calculator
    larflow::reco::NuVertexFlashPrediction predictor;
    larflow::reco::SinkhornFlashDivergence sinkhorn_calc;
    
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
    
    // Process entries
    for (int entry = start_entry; entry < end_entry; entry++) {
        
        if (verbose || entry % 100 == 0) {
            std::cout << "Processing entry " << entry << " / " << end_entry - 1 << std::endl;
        }
        
        // Clear result structure
        memset(&result, 0, sizeof(result));
        result.entry = entry;
        
        // Read data from files
        kps_tree->GetEntry(entry);
        ioman.read_entry(entry);
        ioll.go_to(entry);
        
        // Get ADC images
        auto ev_img = (larcv::EventImage2D*)(ioman.get_data(larcv::kProductImage2D, "wire"));
        if (!ev_img || ev_img->Image2DArray().size() < 3) {
            std::cerr << "Warning: Cannot get ADC images for entry " << entry << std::endl;
            output_tree->Fill();
            continue;
        }
        
        std::vector<larcv::Image2D> adc_v;
        for (size_t p = 0; p < 3; p++) {
            adc_v.push_back(ev_img->Image2DArray()[p]);
        }
        
        // Get observed opflash
        auto ev_opflash = (larlite::event_opflash*)(ioll.get_data(larlite::data::kOpFlash, "simpleFlashBeam"));
        
        result.has_flash = (ev_opflash && ev_opflash->size() > 0);
        
        if (result.has_flash) {
            // Use the first flash (highest PE)
            const auto& flash = ev_opflash->at(0);
            result.obs_total_pe = flash.TotalPE();
            result.obs_time = flash.Time();
            
            for (int pmt = 0; pmt < 32; pmt++) {
                result.obs_pe_per_pmt[pmt] = flash.PE(pmt);
            }
        }
        
        // Process vertex candidates
        result.has_vertex = (nuvetoed_v && nuvetoed_v->size() > 0);
        
        if (!result.has_vertex) {
            // No vertices - fill with defaults and continue
            output_tree->Fill();
            continue;
        }
        
        // Process each vertex candidate
        for (size_t vtx_idx = 0; vtx_idx < nuvetoed_v->size(); vtx_idx++) {
            
            result.vertex_idx = vtx_idx;
            const auto& vertex_candidate = nuvetoed_v->at(vtx_idx);
            
            // Prediction with all particles
            try {
                auto predicted_flash_all = predictor.predictFlash(
                    vertex_candidate,
                    adc_v,
                    adc_threshold,
                    true,   // use_trilinear
                    false   // primary_prongs_only = false (all particles)
                );
                
                result.prediction_success_all = true;
                result.pred_total_pe_all = predictor.getTotalPredictedPE();
                result.n_tracks_all = predictor.getNumTracksProcessed();
                result.n_showers_all = predictor.getNumShowersProcessed();
                result.total_charge_all = predictor.getTotalChargeCollected();
                result.total_photons_all = predictor.getTotalPhotonsEmitted();
                
                // Get per-PMT predictions
                const auto& pe_per_pmt_all = predictor.getPredictedPE();
                for (int pmt = 0; pmt < 32; pmt++) {
                    auto it = pe_per_pmt_all.find(pmt);
                    result.pred_pe_per_pmt_all[pmt] = (it != pe_per_pmt_all.end()) ? it->second : 0.0;
                }
                
            } catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "Warning: Flash prediction (all) failed for entry " << entry 
                              << ", vertex " << vtx_idx << ": " << e.what() << std::endl;
                }
                result.prediction_success_all = false;
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
                
                result.prediction_success_primary = true;
                result.pred_total_pe_primary = predictor.getTotalPredictedPE();
                result.n_tracks_primary = predictor.getNumTracksProcessed();
                result.n_showers_primary = predictor.getNumShowersProcessed();
                result.total_charge_primary = predictor.getTotalChargeCollected();
                result.total_photons_primary = predictor.getTotalPhotonsEmitted();
                
                // Get per-PMT predictions
                const auto& pe_per_pmt_primary = predictor.getPredictedPE();
                for (int pmt = 0; pmt < 32; pmt++) {
                    auto it = pe_per_pmt_primary.find(pmt);
                    result.pred_pe_per_pmt_primary[pmt] = (it != pe_per_pmt_primary.end()) ? it->second : 0.0;
                }
                
            } catch (const std::exception& e) {
                if (verbose) {
                    std::cerr << "Warning: Flash prediction (primary) failed for entry " << entry 
                              << ", vertex " << vtx_idx << ": " << e.what() << std::endl;
                }
                result.prediction_success_primary = false;
            }
            
            // Calculate metrics if we have both prediction and observation
            if (result.has_flash && result.obs_total_pe > 0) {
                
                // PE differences and ratios
                if (result.prediction_success_all) {
                    result.pe_diff_all = result.pred_total_pe_all - result.obs_total_pe;
                    result.pe_ratio_all = result.pred_total_pe_all / result.obs_total_pe;
                    
                    // Calculate Sinkhorn divergences for all particles
                    std::vector<float> pred_pe_vec_all(result.pred_pe_per_pmt_all, result.pred_pe_per_pmt_all + 32);
                    std::vector<float> obs_pe_vec(result.obs_pe_per_pmt, result.obs_pe_per_pmt + 32);
                    
                    for (int i = 0; i < 3; i++) {
                        try {
                            result.sinkhorn_div_all[i] = sinkhorn_calc.calculateDivergence(
                                pred_pe_vec_all,
                                obs_pe_vec,
                                sinkhorn_regularizations[i],
                                100,    // max_iterations
                                1e-6    // tolerance
                            );
                        } catch (const std::exception& e) {
                            if (verbose) {
                                std::cerr << "Warning: Sinkhorn calculation failed (all, reg=" 
                                          << sinkhorn_regularizations[i] << "): " << e.what() << std::endl;
                            }
                            result.sinkhorn_div_all[i] = -1.0; // Invalid value
                        }
                    }
                }
                
                if (result.prediction_success_primary) {
                    result.pe_diff_primary = result.pred_total_pe_primary - result.obs_total_pe;
                    result.pe_ratio_primary = result.pred_total_pe_primary / result.obs_total_pe;
                    
                    // Calculate Sinkhorn divergences for primary particles
                    std::vector<float> pred_pe_vec_primary(result.pred_pe_per_pmt_primary, result.pred_pe_per_pmt_primary + 32);
                    std::vector<float> obs_pe_vec(result.obs_pe_per_pmt, result.obs_pe_per_pmt + 32);
                    
                    for (int i = 0; i < 3; i++) {
                        try {
                            result.sinkhorn_div_primary[i] = sinkhorn_calc.calculateDivergence(
                                pred_pe_vec_primary,
                                obs_pe_vec,
                                sinkhorn_regularizations[i],
                                100,    // max_iterations
                                1e-6    // tolerance
                            );
                        } catch (const std::exception& e) {
                            if (verbose) {
                                std::cerr << "Warning: Sinkhorn calculation failed (primary, reg=" 
                                          << sinkhorn_regularizations[i] << "): " << e.what() << std::endl;
                            }
                            result.sinkhorn_div_primary[i] = -1.0; // Invalid value
                        }
                    }
                }
            }
            
            // Fill output tree for this vertex
            output_tree->Fill();
        }
    }
    
    // Write output and cleanup
    output_tfile->cd();
    output_tree->Write();
    
    std::cout << "\nProcessing complete!" << std::endl;
    std::cout << "Output entries written: " << output_tree->GetEntries() << std::endl;
    
    // Print summary statistics
    output_tree->Draw("pred_total_pe_all", "prediction_success_all", "goff");
    if (output_tree->GetSelectedRows() > 0) {
        double* pred_vals = output_tree->GetV1();
        double sum = 0, min_val = 1e9, max_val = -1e9;
        for (Long64_t i = 0; i < output_tree->GetSelectedRows(); i++) {
            sum += pred_vals[i];
            if (pred_vals[i] < min_val) min_val = pred_vals[i];
            if (pred_vals[i] > max_val) max_val = pred_vals[i];
        }
        std::cout << "\nPredicted PE (all particles) statistics:" << std::endl;
        std::cout << "  Mean: " << sum / output_tree->GetSelectedRows() << std::endl;
        std::cout << "  Min: " << min_val << std::endl;
        std::cout << "  Max: " << max_val << std::endl;
    }
    
    output_tfile->Close();
    reco_tfile->Close();
    ioman.finalize();
    ioll.close();
    
    std::cout << "\nOutput saved to: " << output_file << std::endl;
    
    return 0;
}