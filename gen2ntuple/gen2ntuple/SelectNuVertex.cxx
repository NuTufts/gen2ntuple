#include "SelectNuVertex.h"
#include <iostream>

namespace gen2ntuple {

    int SelectNuVertex::selectVertex(larcv::IOManager* larcv_io, 
                                    larlite::storage_manager* larlite_io, 
                                    gen2ntuple::EventData* event_data,
                                    gen2ntuple::RecoData* reco_data) {
        
        if (!reco_data || !reco_data->nuvtx_v || reco_data->nuvtx_v->empty()) {
            return -1;
        }
        
        const auto& nuvtx_v = *(reco_data->nuvtx_v);
        
        int best_index = -1;
        float best_score = -1.0;
        
        // Loop over all vertex candidates
        for (size_t ivtx = 0; ivtx < nuvtx_v.size(); ivtx++) {
            const auto& vtx = nuvtx_v[ivtx];
            
            // Only consider neutrino keypoints (type == 0)
            if (vtx.keypoint_type != 0) {
                continue;
            }
            
            // Check if this has the highest neutrino score
            if (vtx.netNuScore > best_score) {
                best_score = vtx.netNuScore;
                best_index = static_cast<int>(ivtx);
            }
        }
        
        if (best_index >= 0) {
            std::cout << "SelectNuVertex: Selected vertex " << best_index 
                      << " with neutrino score " << best_score << std::endl;
        }
        
        return best_index;
    }

}