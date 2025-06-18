#include "HighestKPRankWithVisEnergy.h"
#include <iostream>

namespace gen2ntuple {

    int HighestKPRankWithVisEnergy::selectVertex(larcv::IOManager* larcv_io, 
                                                 larlite::storage_manager* larlite_io, 
                                                 gen2ntuple::EventData* event_data,
                                                 gen2ntuple::RecoData* reco_data) {
        
        if (!reco_data || !reco_data->nuvtx_v || reco_data->nuvtx_v->empty()) {
            return -1;
        }
        
        const auto& nuvtx_v = *(reco_data->nuvtx_v);
        const auto& nusel_v = *(reco_data->nusel_v);
        
        // Note: In the Python version, nuselvar_v contains per-vertex selection variables
        // including approx_vis_energy_MeV. We need to calculate or access this.
        // For now, I'll assume this is calculated elsewhere and stored in the vertex candidates
        // or we can use the sum of track and shower energies as a proxy
        
        // First, look for neutrino keypoints (type == 0)
        int best_nu_index = -1;
        float best_nu_vis_energy = -1.0;
        
        for (size_t ivtx = 0; ivtx < nuvtx_v.size(); ivtx++) {
            const auto& vtx = nuvtx_v[ivtx];
            
            // Only consider neutrino keypoints
            if (vtx.keypoint_type != 0) {
                continue;
            }
            
            // count primaries
            int nprim_showers = 0;
            int nprim_tracks  = 0;
            for (int ishower=0; ishower<(int)vtx.shower_v.size(); ishower++) {
                if ( vtx.shower_isSecondary_v.at(ishower)==0 )
                    nprim_showers++;
            }
            for (int itrack=0; itrack<(int)vtx.track_v.size(); itrack++) {
                if ( vtx.track_isSecondary_v.at(itrack)==0 ){
                    nprim_tracks++;
                }
            }

            // Check minimum requirements
            if ( primary_only_ ) {
                if ( nprim_showers<min_num_showers_ )
                    continue;
                if ( nprim_tracks<min_num_tracks_ )
                    continue;
            }
            else {
                if (static_cast<int>(vtx.shower_v.size()) < min_num_showers_) {
                    continue;
                }
                if (static_cast<int>(vtx.track_v.size()) < min_num_tracks_) {
                    continue;
                }
            }
            
            // Get visible energy from NuSelectionVariables
            float vis_energy = 0.0;
            if (ivtx < nusel_v.size()) {
                const auto& selvar = nusel_v[ivtx];
                vis_energy = selvar.approx_vis_energy_MeV;
            }
            
            // Check if this has the highest visible energy among nu vertices
            if (vis_energy > best_nu_vis_energy) {
                best_nu_vis_energy = vis_energy;
                best_nu_index = static_cast<int>(ivtx);
            }
        }
        
        // If we found a neutrino vertex, return it
        if (best_nu_index >= 0) {
            std::cout << "HighestKPRankWithVisEnergy: Selected neutrino vertex " << best_nu_index 
                      << " with visible energy " << best_nu_vis_energy << " MeV" << std::endl;
            return best_nu_index;
        }
        
        // If no neutrino vertex found, look for non-neutrino vertices with highest visible energy
        int best_index = -1;
        float best_vis_energy = -1.0;
        
        for (size_t ivtx = 0; ivtx < nuvtx_v.size(); ivtx++) {
            const auto& vtx = nuvtx_v[ivtx];
            
            // Skip neutrino keypoints this time
            if (vtx.keypoint_type == 0) {
                continue;
            }
            
            // count primaries
            int nprim_showers = 0;
            int nprim_tracks  = 0;
            for (int ishower=0; ishower<(int)vtx.shower_v.size(); ishower++) {
                if ( vtx.shower_isSecondary_v.at(ishower)==0 )
                    nprim_showers++;
            }
            for (int itrack=0; itrack<(int)vtx.track_v.size(); itrack++) {
                if ( vtx.track_isSecondary_v.at(itrack)==0 ){
                    nprim_tracks++;
                }
            }

            // Check minimum requirements
            if ( primary_only_ ) {
                if ( nprim_showers<min_num_showers_ )
                    continue;
                if ( nprim_tracks<min_num_tracks_ )
                    continue;
            }
            else {
                if (static_cast<int>(vtx.shower_v.size()) < min_num_showers_) {
                    continue;
                }
                if (static_cast<int>(vtx.track_v.size()) < min_num_tracks_) {
                    continue;
                }
            }
            
            // Get visible energy from NuSelectionVariables
            float vis_energy = 0.0;
            if (ivtx < nusel_v.size()) {
                const auto& selvar = nusel_v[ivtx];
                vis_energy = selvar.approx_vis_energy_MeV;
            }
            
            // Check if this has the highest visible energy
            if (vis_energy > best_vis_energy) {
                best_vis_energy = vis_energy;
                best_index = static_cast<int>(ivtx);
            }
        }
        
        if (best_index >= 0) {
            std::cout << "HighestKPRankWithVisEnergy: Selected non-neutrino vertex " << best_index 
                      << " with visible energy " << best_vis_energy << " MeV" << std::endl;
        }
        
        return best_index;
    }

}