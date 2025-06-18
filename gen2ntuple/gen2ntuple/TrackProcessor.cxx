#include "TrackProcessor.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// LArLite includes
#include "DataFormat/storage_manager.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcpart.h"

// LArCV includes
#include "larcv/core/DataFormat/IOManager.h"

namespace gen2ntuple {

TrackProcessor::TrackProcessor() 
    : is_mc_(false), vertex_x_(0), vertex_y_(0), vertex_z_(0) {
}

bool TrackProcessor::processEvent(larlite::storage_manager* larlite_io,
                                 larcv::IOManager* larcv_io,
                                 EventData* event_data) {
    
    if (!event_data) {
        std::cerr << "TrackProcessor: EventData pointer is null" << std::endl;
        return false;
    }
    
    // Set vertex position from event data
    setVertexPosition(event_data->vtxX, event_data->vtxY, event_data->vtxZ);
    
    // Extract track information from LArLite
    if (!extractTrackInfo(larlite_io, event_data)) {
        // No tracks found - set defaults
        event_data->nTracks = 0;
        return true;
    }
    
    return true;
}

bool TrackProcessor::extractTrackInfo(larlite::storage_manager* larlite_io,
                                     EventData* event_data) {
    
    // Get track collection from LArLite
    auto ev_track = larlite_io->get_data<larlite::event_track>("nuvertex");
    if (!ev_track || ev_track->size() == 0) {
        event_data->nTracks = 0;
        return true;
    }
    
    int n_tracks = std::min((int)ev_track->size(), EventData::MAX_TRACKS);
    event_data->nTracks = n_tracks;
    
    // Process each track
    for (int i = 0; i < n_tracks; i++) {
        const auto& track = ev_track->at(i);
        
        // Calculate geometric properties
        if (!calculateTrackGeometry(track, i, event_data)) {
            std::cerr << "TrackProcessor: Failed to calculate geometry for track " << i << std::endl;
            continue;
        }
        
        // Calculate angular properties
        if (!calculateTrackAngles(track, i, event_data)) {
            std::cerr << "TrackProcessor: Failed to calculate angles for track " << i << std::endl;
            continue;
        }
        
        // Calculate energy
        if (!calculateTrackEnergy(track, i, event_data)) {
            std::cerr << "TrackProcessor: Failed to calculate energy for track " << i << std::endl;
            continue;
        }
        
        // Set placeholder values for CNN scores (will be filled by CNN integration later)
        event_data->trackElScore[i] = 0.0f;
        event_data->trackPhScore[i] = 0.0f;
        event_data->trackMuScore[i] = 0.5f; // Default to muon-like
        event_data->trackPiScore[i] = 0.0f;
        event_data->trackPrScore[i] = 0.0f;
        event_data->trackPID[i] = PID_MUON; // Default classification
        
        // Set placeholder values for CNN quality metrics
        event_data->trackComp[i] = 0.5f;
        event_data->trackPurity[i] = 0.5f;
        event_data->trackProcess[i] = 0;
        
        // Set placeholder values for origin scores
        event_data->trackPrimaryScore[i] = 0.8f;
        event_data->trackFromNeutralScore[i] = 0.1f;
        event_data->trackFromChargedScore[i] = 0.1f;
        
        // Set basic info
        event_data->trackIsSecondary[i] = isSecondaryTrack(track);
        event_data->trackNHits[i] = track.NumberTrajectoryPoints(); // Approximation
        
        // Set placeholder charge values (would need hit cluster information)
        event_data->trackCharge[i] = 1000.0f; // Placeholder
        event_data->trackChargeFrac[i] = 0.1f;
        event_data->trackHitFrac[i] = 0.1f;
        
        // Truth matching for MC
        if (is_mc_) {
            if (!performTruthMatching(larlite_io, i, event_data)) {
                // Set default truth values
                event_data->trackTruePDG[i] = 0;
                event_data->trackTrueTID[i] = -1;
                event_data->trackTrueMID[i] = -1;
                event_data->trackTrueE[i] = -1.0f;
                event_data->trackTrueComp[i] = 0.0f;
                event_data->trackTruePurity[i] = 0.0f;
            }
        }
    }
    
    return true;
}

bool TrackProcessor::calculateTrackGeometry(const larlite::track& track, int track_idx,
                                           EventData* event_data) {
    
    if (track.NumberTrajectoryPoints() < 2) {
        // Not enough points to calculate geometry
        return false;
    }
    
    // Start position
    auto start_pos = track.Vertex();
    event_data->trackStartPosX[track_idx] = start_pos.X();
    event_data->trackStartPosY[track_idx] = start_pos.Y();
    event_data->trackStartPosZ[track_idx] = start_pos.Z();
    
    // End position
    auto end_pos = track.End();
    event_data->trackEndPosX[track_idx] = end_pos.X();
    event_data->trackEndPosY[track_idx] = end_pos.Y();
    event_data->trackEndPosZ[track_idx] = end_pos.Z();
    
    // Start direction
    auto start_dir = track.VertexDirection();
    event_data->trackStartDirX[track_idx] = start_dir.X();
    event_data->trackStartDirY[track_idx] = start_dir.Y();
    event_data->trackStartDirZ[track_idx] = start_dir.Z();
    
    // Distance to vertex
    event_data->trackDistToVtx[track_idx] = calculateDistanceToVertex(track);
    
    return true;
}

bool TrackProcessor::calculateTrackAngles(const larlite::track& track, int track_idx,
                                         EventData* event_data) {
    
    // Calculate cosine of angle with beam direction (z-axis)
    event_data->trackCosTheta[track_idx] = calculateCosTheta(track);
    
    // Calculate cosine of angle with gravity direction (negative y-axis)
    event_data->trackCosThetaY[track_idx] = calculateCosThetaY(track);
    
    return true;
}

bool TrackProcessor::calculateTrackEnergy(const larlite::track& track, int track_idx,
                                         EventData* event_data) {
    
    float track_length = calculateTrackLength(track);
    
    // For now, use range-based energy for muons (placeholder)
    // In the full implementation, this would use the CNN PID scores
    float muon_energy = calculateMuonEnergy(track);
    float pion_energy = calculatePionRangeEnergy(track_length);
    
    // Default to muon energy (will be refined with CNN)
    event_data->trackRecoE[track_idx] = muon_energy;
    
    return true;
}

bool TrackProcessor::performTruthMatching(larlite::storage_manager* larlite_io,
                                         int track_idx, EventData* event_data) {
    
    // This is a simplified version - the full implementation would use
    // sophisticated pixel-level truth matching from the Python code
    
    // Get MC track information
    auto ev_mctrack = larlite_io->get_data<larlite::event_mctrack>("mcreco");
    if (!ev_mctrack || ev_mctrack->size() == 0) {
        return false;
    }
    
    // Simple truth matching based on closest start position
    // (In practice, would use more sophisticated matching)
    float min_dist = 999999.0f;
    int best_match = -1;
    
    auto reco_start = TVector3(event_data->trackStartPosX[track_idx],
                              event_data->trackStartPosY[track_idx], 
                              event_data->trackStartPosZ[track_idx]);
    
    for (size_t i = 0; i < ev_mctrack->size(); i++) {
        const auto& mctrack = ev_mctrack->at(i);
        if (mctrack.size() == 0) continue;
        
        auto mc_start = mctrack.at(0).Position().Vect();
        float dist = (reco_start - mc_start).Mag();
        
        if (dist < min_dist) {
            min_dist = dist;
            best_match = i;
        }
    }
    
    if (best_match >= 0 && min_dist < 5.0f) { // 5cm matching threshold
        const auto& mctrack = ev_mctrack->at(best_match);
        
        event_data->trackTruePDG[track_idx] = mctrack.PdgCode();
        event_data->trackTrueTID[track_idx] = mctrack.TrackID();
        event_data->trackTrueMID[track_idx] = mctrack.MotherTrackID();
        event_data->trackTrueE[track_idx] = mctrack.at(0).E() * 1000.0f; // Convert to MeV
        
        // Placeholder quality metrics
        event_data->trackTrueComp[track_idx] = 0.8f;
        event_data->trackTruePurity[track_idx] = 0.8f;
        
        return true;
    }
    
    return false;
}

float TrackProcessor::calculateTrackLength(const larlite::track& track) const {
    if (track.NumberTrajectoryPoints() < 2) {
        return 0.0f;
    }
    
    return track.Length();
}

float TrackProcessor::calculateDistanceToVertex(const larlite::track& track) const {
    auto start_pos = track.Vertex();
    
    float dx = start_pos.X() - vertex_x_;
    float dy = start_pos.Y() - vertex_y_;
    float dz = start_pos.Z() - vertex_z_;
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

float TrackProcessor::calculateCosTheta(const larlite::track& track) const {
    auto direction = track.VertexDirection();
    
    // Dot product with beam direction (0, 0, 1)
    float cos_theta = direction.Z() / direction.Mag();
    
    return cos_theta;
}

float TrackProcessor::calculateCosThetaY(const larlite::track& track) const {
    auto direction = track.VertexDirection();
    
    // Dot product with gravity direction (0, -1, 0)
    float cos_theta_y = -direction.Y() / direction.Mag();
    
    return cos_theta_y;
}

std::vector<float> TrackProcessor::getTrackDirection(const larlite::track& track) const {
    auto direction = track.VertexDirection();
    
    return {static_cast<float>(direction.X()), 
            static_cast<float>(direction.Y()), 
            static_cast<float>(direction.Z())};
}

float TrackProcessor::calculateMuonEnergy(const larlite::track& track) const {
    // Placeholder muon energy calculation
    // In practice, this would use sophisticated range-energy relationships
    float length = calculateTrackLength(track);
    
    // Simple range-energy formula for muons (very approximate)
    return length * 2.0f; // MeV per cm (rough approximation)
}

float TrackProcessor::calculatePionRangeEnergy(float track_length) const {
    // Placeholder pion range-energy calculation
    // This would implement the pionRange2T function from the Python code
    
    if (track_length < 1.0f) return 0.0f;
    
    // Simple approximation
    return track_length * 1.5f; // MeV per cm
}

bool TrackProcessor::isSecondaryTrack(const larlite::track& track) const {
    // Simple heuristic: tracks starting far from vertex are likely secondary
    float dist_to_vertex = calculateDistanceToVertex(track);
    
    return dist_to_vertex > 5.0f; // 5cm threshold
}

} // namespace gen2ntuple
