#include "ShowerProcessor.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// LArLite includes
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/shower.h"
#include "larlite/DataFormat/vertex.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mcpart.h"

// LArCV includes
#include "larcv/core/DataFormat/IOManager.h"

namespace gen2ntuple {

ShowerProcessor::ShowerProcessor() 
    : is_mc_(false), vertex_x_(0), vertex_y_(0), vertex_z_(0) {
}

bool ShowerProcessor::processEvent(larlite::storage_manager* larlite_io,
                                  larcv::IOManager* larcv_io,
                                  EventData* event_data) {
    
    if (!event_data) {
        std::cerr << "ShowerProcessor: EventData pointer is null" << std::endl;
        return false;
    }
    
    // Set vertex position from event data
    setVertexPosition(event_data->vtxX, event_data->vtxY, event_data->vtxZ);
    
    // Extract shower information from LArLite
    if (!extractShowerInfo(larlite_io, larcv_io, event_data)) {
        // No showers found - set defaults
        event_data->nShowers = 0;
        return true;
    }
    
    return true;
}

bool ShowerProcessor::extractShowerInfo(larlite::storage_manager* larlite_io,
                                       larcv::IOManager* larcv_io,
                                       EventData* event_data) {
    
    // Get shower collection from LArLite
    auto ev_shower = larlite_io->get_data<larlite::event_shower>("nuvertex");
    if (!ev_shower || ev_shower->size() == 0) {
        event_data->nShowers = 0;
        return true;
    }
    
    int n_showers = std::min((int)ev_shower->size(), EventData::MAX_SHOWERS);
    event_data->nShowers = n_showers;
    
    // Process each shower
    for (int i = 0; i < n_showers; i++) {
        const auto& shower = ev_shower->at(i);
        
        // Calculate geometric properties
        if (!calculateShowerGeometry(shower, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate geometry for shower " << i << std::endl;
            continue;
        }
        
        // Calculate angular properties
        if (!calculateShowerAngles(shower, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate angles for shower " << i << std::endl;
            continue;
        }
        
        // Calculate charge (placeholder - would need cluster information)
        if (!calculateShowerCharge(larlite_io, larcv_io, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate charge for shower " << i << std::endl;
            // Continue with defaults
        }
        
        // Calculate energy
        if (!calculateShowerEnergy(shower, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate energy for shower " << i << std::endl;
            continue;
        }
        
        // Set placeholder values for CNN scores (will be filled by CNN integration later)
        event_data->showerElScore[i] = 0.8f; // Default to electron-like for EM showers
        event_data->showerPhScore[i] = 0.2f;
        event_data->showerMuScore[i] = 0.0f;
        event_data->showerPiScore[i] = 0.0f;
        event_data->showerPrScore[i] = 0.0f;
        event_data->showerPID[i] = PID_ELECTRON; // Default classification
        
        // Set placeholder values for CNN quality metrics
        event_data->showerComp[i] = 0.6f;
        event_data->showerPurity[i] = 0.6f;
        event_data->showerProcess[i] = 0;
        
        // Set placeholder values for origin scores
        event_data->showerPrimaryScore[i] = 0.7f;
        event_data->showerFromNeutralScore[i] = 0.2f;
        event_data->showerFromChargedScore[i] = 0.1f;
        
        // Set basic info
        event_data->showerIsSecondary[i] = isSecondaryShower(shower);
        event_data->showerNHits[i] = static_cast<int>(shower.best_plane()); // Placeholder approximation
        
        // Set placeholder charge values (would need hit cluster information)
        event_data->showerCharge[i] = 2000.0f; // Placeholder - EM showers typically have more charge
        event_data->showerChargeFrac[i] = 0.15f;
        event_data->showerHitFrac[i] = 0.15f;
        
        // Truth matching for MC
        if (is_mc_) {
            if (!performTruthMatching(larlite_io, i, event_data)) {
                // Set default truth values
                event_data->showerTruePDG[i] = 0;
                event_data->showerTrueTID[i] = -1;
                event_data->showerTrueMID[i] = -1;
                event_data->showerTrueE[i] = -1.0f;
                event_data->showerTrueComp[i] = 0.0f;
                event_data->showerTruePurity[i] = 0.0f;
            }
        }
    }
    
    return true;
}

bool ShowerProcessor::calculateShowerGeometry(const larlite::shower& shower, int shower_idx,
                                             EventData* event_data) {
    
    // Start position
    auto start_pos = shower.ShowerStart();
    event_data->showerStartPosX[shower_idx] = start_pos.X();
    event_data->showerStartPosY[shower_idx] = start_pos.Y();
    event_data->showerStartPosZ[shower_idx] = start_pos.Z();
    
    // Start direction
    auto start_dir = shower.Direction();
    event_data->showerStartDirX[shower_idx] = start_dir.X();
    event_data->showerStartDirY[shower_idx] = start_dir.Y();
    event_data->showerStartDirZ[shower_idx] = start_dir.Z();
    
    // Distance to vertex
    event_data->showerDistToVtx[shower_idx] = calculateDistanceToVertex(shower);
    
    return true;
}

bool ShowerProcessor::calculateShowerAngles(const larlite::shower& shower, int shower_idx,
                                           EventData* event_data) {
    
    // Calculate cosine of angle with beam direction (z-axis)
    event_data->showerCosTheta[shower_idx] = calculateCosTheta(shower);
    
    // Calculate cosine of angle with gravity direction (negative y-axis)
    event_data->showerCosThetaY[shower_idx] = calculateCosThetaY(shower);
    
    return true;
}

bool ShowerProcessor::calculateShowerCharge(larlite::storage_manager* larlite_io,
                                           larcv::IOManager* larcv_io,
                                           int shower_idx, EventData* event_data) {
    
    // Placeholder implementation - would need access to associated clusters
    // In the full implementation, this would integrate charge from hits
    // associated with the shower across different planes
    
    // For now, set reasonable defaults
    event_data->showerCharge[shower_idx] = 2000.0f; // Placeholder
    event_data->showerChargeFrac[shower_idx] = 0.15f;
    event_data->showerHitFrac[shower_idx] = 0.15f;
    
    return true;
}

bool ShowerProcessor::calculateShowerEnergy(const larlite::shower& shower, int shower_idx,
                                           EventData* event_data) {
    
    // Calculate different energy estimators
    float electron_energy = calculateElectronEnergy(shower);
    float photon_energy = calculatePhotonEnergy(shower);
    
    // Per-plane energy reconstruction
    float energy_u = calculateShowerEnergyByPlane(shower, 0); // U plane
    float energy_v = calculateShowerEnergyByPlane(shower, 1); // V plane  
    float energy_y = calculateShowerEnergyByPlane(shower, 2); // Y plane
    
    // Store energies
    event_data->showerRecoE[shower_idx] = electron_energy; // Default to electron energy
    event_data->showerRecoEU[shower_idx] = energy_u;
    event_data->showerRecoEV[shower_idx] = energy_v;
    event_data->showerRecoEY[shower_idx] = energy_y;
    
    return true;
}

bool ShowerProcessor::performTruthMatching(larlite::storage_manager* larlite_io,
                                          int shower_idx, EventData* event_data) {
    
    // This is a simplified version - the full implementation would use
    // sophisticated pixel-level truth matching from the Python code
    
    // Get MC shower information
    auto ev_mcshower = larlite_io->get_data<larlite::event_mcshower>("mcreco");
    if (!ev_mcshower || ev_mcshower->size() == 0) {
        return false;
    }
    
    // Simple truth matching based on closest start position
    // (In practice, would use more sophisticated matching)
    float min_dist = 999999.0f;
    int best_match = -1;
    
    auto reco_start = TVector3(event_data->showerStartPosX[shower_idx],
                              event_data->showerStartPosY[shower_idx], 
                              event_data->showerStartPosZ[shower_idx]);
    
    for (size_t i = 0; i < ev_mcshower->size(); i++) {
        const auto& mcshower = ev_mcshower->at(i);
        
        auto mc_start = mcshower.Start().Position().Vect();
        float dist = (reco_start - mc_start).Mag();
        
        if (dist < min_dist) {
            min_dist = dist;
            best_match = i;
        }
    }
    
    if (best_match >= 0 && min_dist < 5.0f) { // 5cm matching threshold
        const auto& mcshower = ev_mcshower->at(best_match);
        
        event_data->showerTruePDG[shower_idx] = mcshower.PdgCode();
        event_data->showerTrueTID[shower_idx] = mcshower.TrackID();
        event_data->showerTrueMID[shower_idx] = mcshower.MotherTrackID();
        event_data->showerTrueE[shower_idx] = mcshower.Start().E() * 1000.0f; // Convert to MeV
        
        // Placeholder quality metrics
        event_data->showerTrueComp[shower_idx] = 0.7f;
        event_data->showerTruePurity[shower_idx] = 0.7f;
        
        return true;
    }
    
    return false;
}

float ShowerProcessor::calculateDistanceToVertex(const larlite::shower& shower) const {
    auto start_pos = shower.ShowerStart();
    
    float dx = start_pos.X() - vertex_x_;
    float dy = start_pos.Y() - vertex_y_;
    float dz = start_pos.Z() - vertex_z_;
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

float ShowerProcessor::calculateCosTheta(const larlite::shower& shower) const {
    auto direction = shower.Direction();
    
    // Dot product with beam direction (0, 0, 1)
    float cos_theta = direction.Z() / direction.Mag();
    
    return cos_theta;
}

float ShowerProcessor::calculateCosThetaY(const larlite::shower& shower) const {
    auto direction = shower.Direction();
    
    // Dot product with gravity direction (0, -1, 0)
    float cos_theta_y = -direction.Y() / direction.Mag();
    
    return cos_theta_y;
}

std::vector<float> ShowerProcessor::getShowerDirection(const larlite::shower& shower) const {
    auto direction = shower.Direction();
    
    return {static_cast<float>(direction.X()), 
            static_cast<float>(direction.Y()), 
            static_cast<float>(direction.Z())};
}

float ShowerProcessor::calculateElectronEnergy(const larlite::shower& shower) const {
    // Get shower energy from LArLite shower object
    // This typically comes from calorimetric reconstruction
    
    // Try to get energy from the shower object using vector method
    auto energy_vec = shower.Energy_v();
    if (energy_vec.size() > 0) {
        // Use the best plane energy (usually the collection plane)
        int best_plane = shower.best_plane();
        if (best_plane >= 0 && best_plane < (int)energy_vec.size()) {
            return energy_vec[best_plane];
        }
        // Fall back to first available energy
        return energy_vec[0];
    }
    
    // Fallback calculation based on dE/dx if no energy available
    auto dedx_vec = shower.dEdx_v();
    if (dedx_vec.size() > 0) {
        // Rough approximation: assume shower length and use dE/dx
        float avg_dedx = 0.0f;
        for (const auto& dedx : dedx_vec) {
            avg_dedx += dedx;
        }
        avg_dedx /= dedx_vec.size();
        
        // Estimate length from shower object (placeholder)
        float estimated_length = 10.0f; // cm, would need proper calculation
        return avg_dedx * estimated_length;
    }
    
    return 100.0f; // Default fallback energy in MeV
}

float ShowerProcessor::calculatePhotonEnergy(const larlite::shower& shower) const {
    // For photons, energy calculation is similar to electrons
    // but may have different calibration factors
    
    float electron_energy = calculateElectronEnergy(shower);
    
    // Apply photon-specific correction factor (if any)
    // This would be calibrated based on MC studies
    float photon_correction = 1.0f; // No correction for now
    
    return electron_energy * photon_correction;
}

float ShowerProcessor::calculateShowerEnergyByPlane(const larlite::shower& shower, int plane) const {
    // Get energy for specific wire plane using vector method
    auto energy_vec = shower.Energy_v();
    
    if (plane >= 0 && plane < (int)energy_vec.size()) {
        return energy_vec[plane];
    }
    
    // Fallback: use dE/dx for this plane if available
    auto dedx_vec = shower.dEdx_v();
    if (plane >= 0 && plane < (int)dedx_vec.size()) {
        float dedx = dedx_vec[plane];
        float estimated_length = 10.0f; // cm, placeholder
        return dedx * estimated_length;
    }
    
    return 50.0f; // Default fallback energy per plane
}

bool ShowerProcessor::isSecondaryShower(const larlite::shower& shower) const {
    // Simple heuristic: showers starting far from vertex are likely secondary
    float dist_to_vertex = calculateDistanceToVertex(shower);
    
    // Also consider energy threshold - very low energy showers might be secondary
    float shower_energy = calculateElectronEnergy(shower);
    
    return (dist_to_vertex > 3.0f) || (shower_energy < 20.0f); // 3cm distance or 20 MeV threshold
}

} // namespace gen2ntuple
