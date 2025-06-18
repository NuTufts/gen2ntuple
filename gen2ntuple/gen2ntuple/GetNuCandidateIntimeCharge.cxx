#include "GetNuCandidateIntimeCharge.h"
#include "larflow/Reco/ClusterImageMask.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include <iostream>
#include <algorithm>
#include <map>

namespace gen2ntuple {

    int GetNuCandidateIntimeCharge::selectVertex(larcv::IOManager* larcv_io, 
                                                 larlite::storage_manager* larlite_io, 
                                                 gen2ntuple::EventData* event_data,
                                                 gen2ntuple::RecoData* reco_data) {
        
        if (!reco_data || !reco_data->nuvtx_v || reco_data->nuvtx_v->empty()) {
            return -1;
        }
        
        const auto& nuvtx_v = *(reco_data->nuvtx_v);
        const auto& nusel_v = *(reco_data->nusel_v);
        
        // Map to store best unreco fraction per keypoint type
        std::map<int, float> kp_intime_reco_frac;
        std::map<int, int> kp_index;
        
        // Initialize maps for keypoint types 0-3
        for (int kp = 0; kp <= 3; kp++) {
            kp_intime_reco_frac[kp] = 10.0; // High initial value
            kp_index[kp] = -1;
        }
        
        float max_intime_reco_frac = 10.0;
        int max_kpindex = -1;
        
        std::cout << "GetNuCandidateIntimeCharge: Processing " << nuvtx_v.size() << " vertices" << std::endl;
        
        // Loop over all vertex candidates
        for (size_t ivtx = 0; ivtx < nuvtx_v.size(); ivtx++) {
            const auto& vtx = nuvtx_v[ivtx];
            
            // Skip vertices with no tracks or showers
            if (vtx.track_v.empty() && vtx.shower_v.empty()) {
                continue;
            }
            
            int kptype = static_cast<int>(vtx.keypoint_type);
            if (kptype > 3) {
                continue;
            }
            
            // Get unreco fraction from NuSelectionVariables
            float unreco_frac = 10.0;  // High default value
            if (ivtx < nusel_v.size()) {
                const auto& selvar = nusel_v[ivtx];
                if (selvar.unreco_fraction_v.size() >= 3) {
                    // Get the middle value (second element) after sorting like in Python
                    std::vector<float> frac_copy = {selvar.unreco_fraction_v[0], 
                                                   selvar.unreco_fraction_v[1], 
                                                   selvar.unreco_fraction_v[2]};
                    std::sort(frac_copy.begin(), frac_copy.end());
                    unreco_frac = frac_copy[1];  // Middle value
                }
            }
            
            // Update best for this keypoint type
            if (unreco_frac < kp_intime_reco_frac[kptype]) {
                kp_intime_reco_frac[kptype] = unreco_frac;
                kp_index[kptype] = static_cast<int>(ivtx);
            }
            
            // Update overall best
            if (unreco_frac < max_intime_reco_frac) {
                max_intime_reco_frac = unreco_frac;
                max_kpindex = static_cast<int>(ivtx);
            }
        }
        
        if (prioritize_by_keypoint_) {
            // Priority order: 0, 3, 1, 2
            std::vector<int> priority_order = {0, 3, 1, 2};
            
            for (int ikp : priority_order) {
                if (kp_index[ikp] != -1) {
                    int selected_idx = kp_index[ikp];
                    float kpscore = nuvtx_v[selected_idx].netScore;
                    std::cout << "GetNuCandidateIntimeCharge: Selected kptype=" << ikp 
                              << " idx=" << selected_idx 
                              << " intime_reco_frac=" << kp_intime_reco_frac[ikp] << std::endl;
                    return selected_idx;
                }
            }
        }
        
        // If not prioritizing by keypoint or no priority vertices found,
        // return the overall best
        if (max_kpindex >= 0) {
            float kpscore = nuvtx_v[max_kpindex].netScore;
            std::cout << "GetNuCandidateIntimeCharge: Selected idx=" << max_kpindex 
                      << " with score=" << kpscore << std::endl;
            return max_kpindex;
        }
        
        return -1;
    }
    
    float GetNuCandidateIntimeCharge::calculateIntimeCharge(const larflow::reco::NuVertexCandidate& vtx,
                                                            larcv::IOManager* larcv_io) const {
        // This is a simplified implementation
        // The full implementation would use ClusterImageMask to mask tracks/showers
        // and calculate the charge not tagged as cosmic
        
        // For now, return a proxy value based on vertex properties
        float charge = 0.0;
        
        // // Sum track lengths as proxy for charge
        // for (const auto& track : vtx.track_v) {
        //     charge += track.pronglen_U + track.pronglen_V + track.pronglen_Y;
        // }
        
        // // Sum shower energies
        // for (const auto& shower : vtx.shower_v) {
        //     charge += shower.trunk_len * 3.0; // Factor for 3 planes
        // }
        
        return charge;
    }

}