#include "FlashPEVertexSelection.h"
#include <iostream>
#include <cmath>
#include <limits>

namespace gen2ntuple {

    int FlashPEVertexSelection::selectVertex(larcv::IOManager* larcv_io, 
                                            larlite::storage_manager* larlite_io, 
                                            gen2ntuple::EventData* event_data,
                                            gen2ntuple::RecoData* reco_data) {
        
        if (!reco_data || !reco_data->nuvtx_v || reco_data->nuvtx_v->empty()) {
            return -1;
        }
        
        const auto& nuvtx_v = *(reco_data->nuvtx_v);
        
        // Note: This method assumes that flash predictions have been calculated
        // for all vertices and stored in the VertexSelector's member variables:
        // _nuvtx_flashpred_v and _nuvtx_sinkhorn_div_v
        // These should be accessible through event_data or passed in some way
        
        int best_index = -1;
        float best_score = std::numeric_limits<float>::max();
        int best_index_passing_cut = -1;
        float best_score_passing_cut = std::numeric_limits<float>::max();
        
        // For now, we'll use a simplified approach
        // In the full implementation, we would access the pre-calculated
        // flash predictions and sinkhorn divergences
        
        std::cout << "FlashPEVertexSelection: Processing " << nuvtx_v.size() << " vertices" << std::endl;
        
        // The observed PE per pmt
        std::vector<float> obs_pe_per_pmt(32,0);
        for (int ipmt=0; ipmt<32; ipmt++) {
            obs_pe_per_pmt[ipmt] = event_data->observedPE[ipmt];
        }
        float observed_pe = event_data->observedPEtotal;

        for (size_t ivtx = 0; ivtx < nuvtx_v.size(); ivtx++) {
            const auto& vtx = nuvtx_v[ivtx];
            
            // Skip vertices with no tracks or showers
            if (vtx.track_v.empty() && vtx.shower_v.empty()) {
                continue;
            }
            
            // Get flash data if available
            float sinkhorn_div = 1.0; // Default value
            float predicted_pe = 0.0;
            
            // retrieve sinkhown divergence and predictions
            if (ivtx < sinkhorn_divs_.size()) {
                sinkhorn_div = sinkhorn_divs_[ivtx];
            }
            
            std::vector<float> predicted_pe_per_pmt;
            if (ivtx < flash_predictions_.size()) {
                // Sum predicted PE across all PMTs
                predicted_pe_per_pmt = flash_predictions_.at(ivtx);
                for (float pe : flash_predictions_[ivtx]) {
                    predicted_pe += pe;
                }
            }

            // fractional error
            float frac_pe_error = (predicted_pe - observed_pe) / (0.1+observed_pe);
            
            // Calculate the combined score
            float score = calculateFlashScore(sinkhorn_div, predicted_pe, observed_pe);
            
            if (score < best_score) {
                best_score = score;
                best_index = static_cast<int>(ivtx);
            }

            float m = -(70.0/3);
            float z = sinkhorn_div - m*(frac_pe_error+1.0) + 70.0;

            if ( z<=0.0 && score < best_score_passing_cut ) {
                best_score_passing_cut = score;
                best_index_passing_cut = static_cast<int>(ivtx);
            }

        }
        
        // use more restrictive selection, if possible
        if (best_index_passing_cut >= 0) {
            best_index = best_index_passing_cut;
        }

        if (best_index>=0)  {
            std::cout << "FlashPEVertexSelection: Selected vertex " << best_index 
                      << " with flash match score " << best_score << std::endl;
        }
        
        return best_index;
    }
    
    float FlashPEVertexSelection::calculateFlashScore(float sinkhorn_div, 
                                                     float pred_pe, 
                                                     float obs_pe) const {
        // Calculate fractional PE error
        float frac_pe_error = (pred_pe - obs_pe) / (0.1+obs_pe);
        
        // Combined score: sinkhorn_divergence * |1 + fractional_pe_error|
        float score = sinkhorn_div * std::abs(1.0 + frac_pe_error);
        
        return score;
    }

}