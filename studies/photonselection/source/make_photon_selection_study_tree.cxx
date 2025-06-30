#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// larcv includes
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// larlite includes
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/opflash.h"
#include "larlite/DataFormat/mctruth.h"

// larflow includes
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/NuVertexFlashPrediction.h"
#include "larflow/Reco/SinkhornFlashDivergence.h"

// ublarcvapp includes
#include "ublarcvapp/MCTools/MCPixelPGraph.h"

// Function declarations
float calculateSinkhornDivergence(const std::vector<float>& pred_pe, 
                                  const std::vector<float>& obs_pe);
float calculatePhotonCharge(const larflow::reco::NuVertexCandidate& vtx);
float calculatePhotonNumHits(const larflow::reco::NuVertexCandidate& vtx);
void calculateTotalPixelSum(const larflow::reco::NuVertexCandidate& vtx, 
                           float totalpixelsum[3]);
float calculateDistance(const std::vector<float>& pos1, const std::vector<float>& pos2);
bool isPrimaryPhoton(const ublarcvapp::mctools::MCPixelPGraph::Node_t& node, 
                     const std::vector<float>& true_nu_vtx);
std::vector<float> getTrueNuVertex(larlite::storage_manager& larlite_mgr);
float getPhotonEdepDistance(ublarcvapp::mctools::MCPixelPGraph& mcpg,
                           const std::vector<float>& reco_vtx,
                           std::vector<float>& edep_pos);
std::vector<float> calculateProngPixelSum( larcv::IOManager& larcv_io,
                                           const larflow::reco::NuVertexCandidate& nuvtx,
					   std::vector<float>& leading_shower_pixsum_v,
                                           bool primary_only );
float dwall( float x, float y, float z );

int main( int nargs, char** argv )
{
    std::cout << "Make Photon Selection Study Tree" << std::endl;
    
    // Parse arguments
    if (nargs < 4) {
        std::cout << "Usage: " << argv[0] << " <dlmerged_file> <reco_file> <output_file> [--tickbackward]" << std::endl;
        std::cout << "  dlmerged_file: Input dlmerged ROOT file" << std::endl;
        std::cout << "  reco_file: Input reco ROOT file with KPSRecoManagerTree" << std::endl;
        std::cout << "  output_file: Output ROOT file for study tree" << std::endl;
        std::cout << "  --tickbackward: Optional flag for TickBackward mode" << std::endl;
        std::cout << "  --printmcpg: Optional flag to print MCPixelPGraph info for simulated particles in each event." << std::endl;
        return 1;
    }
    
    std::string dlmerged_file = argv[1];
    std::string reco_file = argv[2];
    std::string output_file = argv[3];
    bool tickbackward_mode = false;
    bool print_mcpg_info = false;
    
    // Check for optional tickbackward flag
    for (int i = 4; i < nargs; i++) {
        if (std::string(argv[i]) == "--tickbackward") {
            tickbackward_mode = true;
            std::cout << "TickBackward mode enabled" << std::endl;
        }
        else if ( std::string(argv[i])== "--printmcpg" ) {
            print_mcpg_info =  true;
            std::cout << "Print MCPixelPGraph information for each entry" << std::endl;
        }

    }
    
    std::cout << "Input dlmerged file: " << dlmerged_file << std::endl;
    std::cout << "Input reco file: " << reco_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    
    // Open files
    std::cout << "Opening files..." << std::endl;
    
    // Open dlmerged file with larcv IOManager for image-like data
    larcv::IOManager larcv_ioman(larcv::IOManager::kREAD, "IOManager", larcv::IOManager::kTickBackward);
    larcv_ioman.add_in_file(dlmerged_file);
    larcv_ioman.reverse_all_products();
    larcv_ioman.initialize();
    
    // Open dlmerged file with larlite storage_manager for particle data
    larlite::storage_manager larlite_mgr(larlite::storage_manager::kREAD);
    larlite_mgr.add_in_filename(dlmerged_file);
    larlite_mgr.open();
    
    // Open reco file with ROOT TFile
    TFile* reco_tfile = TFile::Open(reco_file.c_str(), "READ");
    if (!reco_tfile || reco_tfile->IsZombie()) {
        std::cerr << "Error: Could not open reco file: " << reco_file << std::endl;
        return 1;
    }
    
    // Setup KPSRecoManagerTree
    std::cout << "Setting up KPSRecoManagerTree..." << std::endl;
    TTree* kps_tree = (TTree*)reco_tfile->Get("KPSRecoManagerTree");
    if (!kps_tree) {
        std::cerr << "Error: Could not find KPSRecoManagerTree in reco file" << std::endl;
        return 1;
    }
    
    // Setup branches for KPSRecoManagerTree
    std::vector<larflow::reco::NuVertexCandidate>* nu_vetoed_v = nullptr;
    std::vector<larflow::reco::NuSelectionVariables>* nu_sel_v = nullptr;
    int run, subrun, event;
    
    kps_tree->SetBranchAddress("nuvetoed_v", &nu_vetoed_v);
    kps_tree->SetBranchAddress("nu_sel_v", &nu_sel_v);
    kps_tree->SetBranchAddress("run", &run);
    kps_tree->SetBranchAddress("subrun", &subrun);
    kps_tree->SetBranchAddress("event", &event);
    
    // Get number of entries
    int nentries_larcv = larcv_ioman.get_n_entries();
    int nentries_kps = kps_tree->GetEntries();
    
    std::cout << "Number of entries in larcv: " << nentries_larcv << std::endl;
    std::cout << "Number of entries in KPS tree: " << nentries_kps << std::endl;
    
    if (nentries_larcv != nentries_kps) {
        std::cerr << "Warning: Number of entries mismatch between files" << std::endl;
        return 1;
    }
    
    int nentries = std::min(nentries_larcv, nentries_kps);
    
    // Create output TTree
    std::cout << "Creating output tree..." << std::endl;
    TFile* output_tfile = TFile::Open(output_file.c_str(), "RECREATE");
    TTree* output_tree = new TTree("PhotonSelectionTree", "Photon selection study data");
    
    // Output tree variables
    int out_run, out_subrun, out_event, out_entry, out_vertexindex;
    float out_sinkhorndiv, out_totpefracerr, out_dist2truenuvtx, out_dist2photonedep;
    float out_true_edeppixsum[3];
    float out_true_pixelsum[3];
    float out_reco_pixelsum[3];    
    float out_allprong_pixelsum[3];
    float out_true_median_pixsum, out_true_median_edep;
    float out_reco_median_pixsum, out_allprong_median_pixsum;
    float out_photoncharge, out_photonnumhits;
    float out_observed_totpe, out_predicted_totpe;
    float out_dwall_true_nuvtx, out_dwall_true_edep, out_dwall_reco_nuvtx;
    
    // Create branches
    output_tree->Branch("run", &out_run);
    output_tree->Branch("subrun", &out_subrun);
    output_tree->Branch("event", &out_event);
    output_tree->Branch("entry", &out_entry);
    output_tree->Branch("vertexindex",     &out_vertexindex);
    output_tree->Branch("observed_totpe",  &out_observed_totpe );
    output_tree->Branch("predicted_totpe", &out_predicted_totpe );
    output_tree->Branch("sinkhorndiv",     &out_sinkhorndiv);
    output_tree->Branch("totpefracerr",    &out_totpefracerr);
    output_tree->Branch("dist2truenuvtx",  &out_dist2truenuvtx);
    output_tree->Branch("dist2photonedep", &out_dist2photonedep);
    output_tree->Branch("true_pixelsum",   out_true_pixelsum,   "true_pixselsum[3]/F");
    output_tree->Branch("true_edeppixsum", out_true_edeppixsum, "true_edeppixsum[3]/F");
    output_tree->Branch("true_median_pixsum", &out_true_median_pixsum );
    output_tree->Branch("true_median_edep",   &out_true_median_edep );
    output_tree->Branch("reco_pixelsum",       out_reco_pixelsum,  "reco_pixelsum[3]/F");
    output_tree->Branch("reco_median_pixsum",     &out_reco_median_pixsum );
    output_tree->Branch("allprong_pixelsum",      out_allprong_pixelsum,  "allprong_pixelsum[3]/F");
    output_tree->Branch("allprong_median_pixsum", &out_allprong_median_pixsum );
    output_tree->Branch("photoncharge",  &out_photoncharge);
    output_tree->Branch("photonnumhits", &out_photonnumhits);
    output_tree->Branch("dwall_true_nuvtx", &out_dwall_true_nuvtx);
    output_tree->Branch("dwall_true_edep",  &out_dwall_true_edep);
    output_tree->Branch("dwall_reco_nuvtx", &out_dwall_reco_nuvtx);
    
    // Setup flash predictor
    larflow::reco::NuVertexFlashPrediction flash_predictor;
    flash_predictor.set_verbosity( larcv::msg::kINFO );
    flash_predictor.setChargeToPhotonParams(100.0, 23.6e-3, 24000.0, 0.7);
    flash_predictor.setTrackConversionParams(3, 3, 0.3, 0.5);
    flash_predictor.setShowerConversionParams(3, 3);
    float adc_threshold = 10.0;
    
    // Start event loop
    std::cout << "Starting event loop over " << nentries << " entries..." << std::endl;
    //nentries = 10; // remove this -- for debug
    
    for (int ientry = 0; ientry < nentries; ientry++) {
        if (ientry % 100 == 0) {
            std::cout << "Processing entry " << ientry << "/" << nentries << std::endl;
        }
        
        // Load entries
        larcv_ioman.read_entry(ientry);
        larlite_mgr.go_to(ientry);
        kps_tree->GetEntry(ientry);
        
        // Set output event info
        out_run = run;
        out_subrun = subrun;
        out_event = event;
        out_entry = ientry;
        
        // Create MCPixelPGraph for MC truth information
        ublarcvapp::mctools::MCPixelPGraph mcpg;
        try {
            mcpg.buildgraph(larcv_ioman, larlite_mgr);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not build MCPixelPGraph for entry " << ientry << ": " << e.what() << std::endl;
            continue;
        }

        if ( print_mcpg_info )
          mcpg.printGraph( 0, false );
        
        // Get true neutrino vertex
        std::vector<float> true_nu_vtx = getTrueNuVertex(larlite_mgr);
        out_dwall_true_nuvtx = dwall( true_nu_vtx[0], true_nu_vtx[1], true_nu_vtx[2] );
        
        // Check if this event has target neutrino interaction (primary photons)
        bool has_primary_photon = false;
	std::vector<float> photon_pixelsum(3,-1.0); // pixel sum for all of the photon
	float photon_median_pixelsum = -1;
	std::vector<float> edep_pixelsum(3,-1.0);   // pixel sum for the first cluster
	float edep_median_pixelsum = -1;

	int max_nodeidx = -1;
	float max_median_pixelsum = -1;
	int nnodes = mcpg.node_v.size();
        for (int inode=0; inode<nnodes; inode++) {
	  const auto& node = mcpg.node_v.at(inode);
	  if (isPrimaryPhoton(node, true_nu_vtx)) {
	    has_primary_photon = true;
	    std::vector<float> planepixsums = node.pixsum_v;
	    std::cout << "primary photon: trackid=" << node.tid << " E=" << node.E_MeV << " "
		      << "pixsum_v=(" << planepixsums[0] << "," << planepixsums[1] << "," << planepixsums[2] << ")"
		      << std::endl;	    
	    if ( planepixsums.size()>=3 ) {
	      // has pixel sum. sort to find median
	      std::sort( planepixsums.begin(), planepixsums.end() );
	      // check if this median pixelsum is the largest we've seen
	      if ( max_median_pixelsum<0 || max_median_pixelsum<planepixsums[1] ) {
		max_median_pixelsum = planepixsums[1];
		max_nodeidx = inode;
		photon_pixelsum = node.pixsum_v; // store unsorted values
		std::vector<float> plane_edep = mcpg.getTruePhotonTrunkPlanePixelSums( node.tid );
		edep_pixelsum = plane_edep; // copy over unsorted values
		std::sort( plane_edep.begin(), plane_edep.end() ); // sort this vector
		edep_median_pixelsum = plane_edep[1]; // use sorted values to get median edep
	      }
	    }
	  }
        }

	// Pass true photon info to output tree variables
	out_true_median_pixsum = max_median_pixelsum;
	out_true_median_edep   = edep_median_pixelsum;
	for (int p=0; p<3; p++) {
	  out_true_pixelsum[p]   = photon_pixelsum[p];
	  out_true_edeppixsum[p] = edep_pixelsum[p];
	}

        
        // Get observed opflash
        float obs_total_pe = 0.0;
        std::vector<float> obs_pe_per_pmt(32, 0.0);
        try {
            auto ev_opflash = (larlite::event_opflash*)(larlite_mgr.get_data(larlite::data::kOpFlash, "simpleFlashBeam"));
            if (ev_opflash && ev_opflash->size() > 0) {
                auto& flash = ev_opflash->at(0);
                
                int nOpDets = std::min(static_cast<size_t>(32), flash.nOpDets());
                for (int ipmt = 0; ipmt < nOpDets; ipmt++) {
                    obs_pe_per_pmt[ipmt] = flash.PE(ipmt);
                    obs_total_pe += obs_pe_per_pmt[ipmt];
                }

            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not get opflash for entry " << ientry << ": " << e.what() << std::endl;
        }
        out_observed_totpe = obs_total_pe;
        
        // Get ADC images for flash prediction
        auto ev_adc = (larcv::EventImage2D*)larcv_ioman.get_data(larcv::kProductImage2D, "wire");
        const std::vector<larcv::Image2D>& adc_v = ev_adc->as_vector();
        
        // Process each vertex candidate
        if (nu_vetoed_v) {
            for (size_t ivtx = 0; ivtx < nu_vetoed_v->size(); ivtx++) {
                larflow::reco::NuVertexCandidate& vtx = nu_vetoed_v->at(ivtx);
                out_vertexindex = (int)ivtx;
                
                // Calculate flash prediction variables
                try {
                    auto predicted_flash = flash_predictor.predictFlash(
                        vtx,
                        adc_v,
                        adc_threshold,
                        true,  // use trilinear interpolation
                        false  // include all particles
                    );
                    
                    float pred_total_pe = flash_predictor.getTotalPredictedPE();
                    auto pe_per_pmt_map = flash_predictor.getPredictedPE();
                    
                    // Convert map to vector for Sinkhorn calculation
                    std::vector<float> pred_pe_per_pmt(32, 0.0);
                    out_predicted_totpe = 0.;
                    for (int ipmt = 0; ipmt < 32; ipmt++) {
                        if (pe_per_pmt_map.find(ipmt) != pe_per_pmt_map.end()) {
                            pred_pe_per_pmt[ipmt] = pe_per_pmt_map[ipmt];
                            out_predicted_totpe += pred_pe_per_pmt[ipmt];
                        }
                    }
                    
                    // Calculate flash prediction metrics
                    out_sinkhorndiv = calculateSinkhornDivergence(pred_pe_per_pmt, obs_pe_per_pmt);
                    out_totpefracerr = (pred_total_pe - obs_total_pe) / (0.1 + obs_total_pe);
                    
                } catch (const std::exception& e) {
                    std::cerr << "Warning: Flash prediction failed for entry " << ientry << ", vertex " << ivtx << ": " << e.what() << std::endl;
                    out_sinkhorndiv = -999.0;
                    out_totpefracerr = -999.0;
                }
                
                
                std::vector<float> reco_vtx = {vtx.pos[0], vtx.pos[1], vtx.pos[2]};
                out_dwall_reco_nuvtx = dwall( reco_vtx[0], reco_vtx[1], reco_vtx[2] );

                // Calculate distance to true neutrino vertex
                out_dist2truenuvtx = calculateDistance(reco_vtx, true_nu_vtx);
                
                // Calculate distance to photon energy deposition
                std::vector<float> edep_pos(3,0);
                out_dist2photonedep = getPhotonEdepDistance(mcpg, reco_vtx, edep_pos);
                if ( out_dist2photonedep>=0 ) {
                    out_dwall_true_edep = dwall( edep_pos[0], edep_pos[1], edep_pos[2] );
                }
                else {
                    out_dwall_true_edep = -999.0;
                }
                
                // Calculate total pixel sum across planes
		std::vector<float> leadingshower_v(3,-1.0);
                std::vector<float> totalpixelsum = calculateProngPixelSum( larcv_ioman, vtx, leadingshower_v, true );
                for (int p=0; p<3; p++) {
		    out_allprong_pixelsum[p] = totalpixelsum[p];
		    out_reco_pixelsum[p]     = leadingshower_v[p];
		}
                std::sort( totalpixelsum.begin(), totalpixelsum.end() );
                out_allprong_median_pixsum = totalpixelsum[1];
		std::sort( leadingshower_v.begin(), leadingshower_v.end() );
		out_reco_median_pixsum = leadingshower_v[1];
                
                // Calculate photon-related variables
                out_photoncharge  = calculatePhotonCharge(vtx);
                out_photonnumhits = calculatePhotonNumHits(vtx);
                
                // Fill output tree
                output_tree->Fill();
            }
        }
    }
    
    // Write output and cleanup
    std::cout << "Writing output file..." << std::endl;
    output_tfile->cd();
    output_tree->Write();
    output_tfile->Close();
    
    // Cleanup
    reco_tfile->Close();
    larlite_mgr.close();
    larcv_ioman.finalize();
    
    std::cout << "Processing complete! Output written to: " << output_file << std::endl;
    std::cout << "Total entries processed: " << nentries << std::endl;
    
    return 0;
}

// Function implementations
float calculateSinkhornDivergence(const std::vector<float>& pred_pe, 
                                  const std::vector<float>& obs_pe) {
    try {
        larflow::reco::SinkhornFlashDivergence sinkhorn_calc;
        return sinkhorn_calc.calculateDivergence(pred_pe, obs_pe, 1.0, 100, 1e-6);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Sinkhorn calculation failed: " << e.what() << std::endl;
        return -999.0;
    }
}

float calculatePhotonCharge(const larflow::reco::NuVertexCandidate& vtx) {
    float max_shower_charge = -1.0;
    
    for (const auto& shower : vtx.shower_v) {
        float shower_charge = 0.0;
        for (const auto& hit : shower) {
            // Use Y-plane charge (collection plane) at index 25 if available
            if (hit.size() > 25) {
                shower_charge += hit[25];  // Y plane charge
            }
        }
        if (shower_charge > max_shower_charge) {
            max_shower_charge = shower_charge;
        }
    }
    
    return max_shower_charge;
}

float calculatePhotonNumHits(const larflow::reco::NuVertexCandidate& vtx) {
    int max_shower_hits = -1;
    
    for (const auto& shower : vtx.shower_v) {
        int nhits = (int)shower.size();
        if (nhits > max_shower_hits) {
            max_shower_hits = nhits;
        }
    }
    
    return static_cast<float>(max_shower_hits);
}

void calculateTotalPixelSum(const larflow::reco::NuVertexCandidate& vtx, 
                           float totalpixelsum[3]) {
    // Initialize
    for (int i = 0; i < 3; i++) totalpixelsum[i] = 0.0;
    
    // Sum over all tracks - sum charge from track hits if available
    for (const auto& track : vtx.track_v) {
        // For now use trajectory points as proxy since track charge access is complex
        size_t npoints = track.NumberTrajectoryPoints();
        for (int plane = 0; plane < 3; plane++) {
            totalpixelsum[plane] += static_cast<float>(npoints);
        }
    }
    
    // Sum over all showers - sum actual charge from hits
    for (const auto& shower : vtx.shower_v) {
        for (const auto& hit : shower) {
            if (hit.size() > 25) {
                totalpixelsum[0] += hit[23];  // U plane charge
                totalpixelsum[1] += hit[24];  // V plane charge  
                totalpixelsum[2] += hit[25];  // Y plane charge
            }
        }
    }
}

float calculateDistance(const std::vector<float>& pos1, const std::vector<float>& pos2) {
    if (pos1.size() < 3 || pos2.size() < 3) return -999.0;
    
    float dx = pos1[0] - pos2[0];
    float dy = pos1[1] - pos2[1];  
    float dz = pos1[2] - pos2[2];
    
    return sqrt(dx*dx + dy*dy + dz*dz);
}

bool isPrimaryPhoton(const ublarcvapp::mctools::MCPixelPGraph::Node_t& node, 
                     const std::vector<float>& true_nu_vtx) {
    if (node.pid != 22) return false;  // Must be photon
    if (node.origin != 1) return false;  // Must be from neutrino
    if (node.first_tpc_pos.size() < 3) return false;
    if (true_nu_vtx.size() < 3) return false;
    
    // Check distance from neutrino vertex
    float dist = sqrt(pow(node.start[0] - true_nu_vtx[0], 2) +
                      pow(node.start[1] - true_nu_vtx[1], 2) +
                      pow(node.start[2] - true_nu_vtx[2], 2));
    
    return dist < 5.0;  // Within 5 cm of vertex
}

std::vector<float> getTrueNuVertex(larlite::storage_manager& larlite_mgr) {
    std::vector<float> vtx_pos(3, -999.0);
    
    try {
        auto ev_mctruth = (larlite::event_mctruth*)(larlite_mgr.get_data(larlite::data::kMCTruth, "generator"));
        if (ev_mctruth && ev_mctruth->size() > 0) {
            auto& mctruth = ev_mctruth->at(0);
            auto& nu = mctruth.GetNeutrino();
            auto& trajectory = nu.Nu().Trajectory();
            if (trajectory.size() > 0) {
                const auto& first_step = trajectory[0];
                vtx_pos[0] = first_step.X();
                vtx_pos[1] = first_step.Y();
                vtx_pos[2] = first_step.Z();
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not get MC truth: " << e.what() << std::endl;
    }
    
    return vtx_pos;
}

float getPhotonEdepDistance(ublarcvapp::mctools::MCPixelPGraph& mcpg,
                           const std::vector<float>& reco_vtx,
                           std::vector<float>& edep_pos ) {

    float min_distance = -1.0;
    edep_pos.resize(3,0.0);
    
    for (const auto& node : mcpg.node_v) {
        if (node.pid == 22 && node.origin == 1) {  // Primary photon
            try {
                const std::vector<std::vector<float>>& photon_points = 
                    mcpg.getTruePhotonTrunk3DPoints(node.tid);
                    
                for (const auto& point : photon_points) {
                    if (point.size() >= 3 && reco_vtx.size() >= 3) {
                        float dist = sqrt(pow(point[0] - reco_vtx[0], 2) +
                                        pow(point[1] - reco_vtx[1], 2) +
                                        pow(point[2] - reco_vtx[2], 2));
                        if (min_distance < 0 || dist < min_distance) {
                            min_distance = dist;
                            edep_pos[0] = point[0];
                            edep_pos[1] = point[1];
                            edep_pos[2] = point[2];
                        }
                    }
                }
            } catch (const std::exception& e) {
                // Continue if this photon doesn't have trunk points
                continue;
            }
        }
    }
    
    return min_distance;
}

std::vector<float> calculateProngPixelSum( larcv::IOManager& larcv_io,
                                           const larflow::reco::NuVertexCandidate& nuvtx,
					   std::vector<float>& leading_shower_pixsum_v,
                                           bool primary_only ) 
{

    int ntracks  = nuvtx.track_v.size();
    int nshowers = nuvtx.shower_v.size();

    std::vector<float> total_charge(3,0.0);
    leading_shower_pixsum_v.resize(3,-1.0);
    for (int p=0; p<3; p++)
      leading_shower_pixsum_v[p] = -1.0;

    auto ev_img = (larcv::EventImage2D*)larcv_io.get_data("image2d", "wire");
    auto image2Dvec = ev_img->as_vector();
    float threshold = 10.0f; // ADC threshold

    for (int p = 0; p < 3; ++p ) {

        auto const& img = image2Dvec.at(p);
        std::set< std::array<int,3> > pixels_visited; // prevent double counting pixels

        for (int track_idx=0; track_idx<ntracks; track_idx++) {

            float clusterCharge = 0.0f;

            if ( primary_only && nuvtx.track_isSecondary_v[track_idx]==1 ) {
                // is secondary and we are only including primaries
                continue;
            }

            auto const& cluster = nuvtx.track_hitcluster_v.at(track_idx);

            for (auto const& hit : cluster) {

                int row = int((hit.tick - 2400) / 6);
                int col = int(hit.targetwire[p]);
                
                if (row >= 0 && row < (int)img.meta().rows() &&
                    col >= 0 && col < (int)img.meta().cols()) {
                    float pixVal = img.pixel(row, col);
                    if (pixVal >= threshold) {
                        std::array<int,3> pixindex = {p,row,col};
                        if ( pixels_visited.find(pixindex)==pixels_visited.end()) {
                            // not in set
                            clusterCharge += pixVal;
                            pixels_visited.insert( pixindex );
                        }
                    }
                }
            }

            total_charge[p] += clusterCharge;

        }//end of track loop

        
        // loop over showers and get pixelsums for the plane
	int leading_shower_index = -1;
	float leading_shower_pixsum = -1.0;
        for (int shower_idx=0; shower_idx<nshowers; shower_idx++ ) {

            float clusterCharge = 0.0f;

            if ( primary_only && nuvtx.shower_isSecondary_v[shower_idx]==1 ) {
                continue;
            }

            auto const& cluster = nuvtx.shower_v.at(shower_idx);

            for (auto const& hit : cluster) {
                int row = int((hit.tick - 2400) / 6);
                int col = int(hit.targetwire[p]);
                
                if (row >= 0 && row < (int)img.meta().rows() &&
                    col >= 0 && col < (int)img.meta().cols()) {
    
                    float pixVal = img.pixel(row, col);
                    if (pixVal >= threshold) {
                        std::array<int,3> pixindex = {p,row,col};
    
                        if ( pixels_visited.find(pixindex)==pixels_visited.end()) {
                            // not in set
                            clusterCharge += pixVal;
                            pixels_visited.insert( pixindex );
                        }
                    }
                }
            }// loop over hits

	    if ( leading_shower_pixsum<0 || clusterCharge>leading_shower_pixsum ) {
	      leading_shower_pixsum = clusterCharge;
	      leading_shower_index = shower_idx;
	    }

            total_charge[p] += clusterCharge;

        }// loop over showers

	leading_shower_pixsum_v[p] = leading_shower_pixsum;
	
    }// loop over planes

    return total_charge;
}

float dwall( float x, float y, float z ) {
    float xmin = x;
    float xmax = 256.0-x;
    float ymin = y + 116.5; // y - -116.5
    float ymax = 116.5-y;
    float zmin = z;
    float zmax = 1036.0-z;

    std::array<float,6> dists = { xmin, xmax, ymin, ymax, zmin, zmax };

    if ( xmin>=0 && xmax>=0 && ymin>=0 && ymax>=0 && zmin>=0 && zmax>=0 ) {
        // is inside, find min
        float mindist = 1.0e9;
        for (auto& dist : dists ) {
            if ( dist < mindist )
                mindist = dist;
        }
        return mindist;
    }
    else {
        // outside: only negative values
        float mindist = 1.0e9;
        for ( auto& dist : dists ) {
            if ( dist>=0 ) continue; // must be on the "wrong" side of the boundary
            if ( fabs(dist)<mindist )
                mindist = fabs(dist);
        }
        // return a negative value to indicate we are outside the TPC
        return -mindist;

    }
    
    // should never get here
    return 1.0e9;
}
