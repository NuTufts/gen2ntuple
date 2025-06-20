#include "ShowerProcessor.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>

// LArLite includes
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/shower.h"
#include "larlite/DataFormat/track.h"
#include "larlite/DataFormat/vertex.h"
#include "larlite/DataFormat/mcshower.h"
#include "larlite/DataFormat/mcpart.h"
#include "larlite/DataFormat/larflowcluster.h"

// LArFlow includes  
#include "larflow/Reco/NuVertexCandidate.h"

// LArCV includes
#include "larcv/core/DataFormat/IOManager.h"
#include "larcv/core/DataFormat/EventImage2D.h"

// LArPID includes
#include "larpid/interface/LArPIDInterface.h"
#include "larpid/model/TorchModel.h"

// UBDLLee includes
#include "ublarcvapp/ubdllee/dwall.h"

// Gen2Ntuple includes
#include "WCFiducial.h"

namespace gen2ntuple {

ShowerProcessor::ShowerProcessor() 
    : is_mc_(false), vertex_x_(0), vertex_y_(0), vertex_z_(0), larpid_cnn_(nullptr) {
}

bool ShowerProcessor::processEvent(larlite::storage_manager* larlite_io,
                                  larcv::IOManager* larcv_io,
                                  EventData* event_data,
                                  RecoData* reco_data) {
    
    if (!event_data) {
        std::cerr << "ShowerProcessor: EventData pointer is null" << std::endl;
        return false;
    }

    if (!reco_data) {
        std::cerr << "ShowerProcessor: RecoData pointer is null" << std::endl;
        return false;
    }
    
    // Set vertex position from event data
    setVertexPosition(event_data->vtxX, event_data->vtxY, event_data->vtxZ);
    
    // Extract shower information from RecoData
    if (!extractShowerInfo(larlite_io, larcv_io, event_data, reco_data)) {
        // No showers found - set defaults
        event_data->nShowers = 0;
        return true;
    }
    
    return true;
}

bool ShowerProcessor::extractShowerInfo(larlite::storage_manager* larlite_io,
                                       larcv::IOManager* larcv_io,
                                       EventData* event_data,
                                       RecoData* reco_data) {
    
    // Get shower collection from RecoData vertex candidates (similar to tracks)
    int vtxIdx = event_data->vtxIndex;
    
    if (reco_data->nuvtx_v->size() == 0) {
        // no reco vertices in this event. just return.
        std::cout << "ShowerProcessor::extractShowerInfo - no vertex in event" << std::endl;
        event_data->nShowers = 0;
        return true;
    }

    if (vtxIdx < 0 || vtxIdx >= (int)reco_data->nuvtx_v->size()) {
        throw std::runtime_error("ShowerProcessor::extractShowerInfo -- vtxIndex unset -- needs to be selected in VertexSelector.");
    }
    
    auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
    
    int n_showers = std::min((int)nuvtx.shower_v.size(), EventData::MAX_SHOWERS);
    event_data->nShowers = n_showers;
    
    std::cout << "ShowerProcessor::extractShowerInfo - number of Showers = " << n_showers << std::endl;
    
    // Process each shower
    for (int i = 0; i < n_showers; i++) {
        const auto& showercluster = nuvtx.shower_v.at(i);  // shower_v contains larflowcluster objects
        const auto& shower_trunk = nuvtx.shower_trunk_v.at(i);  // trunk for direction/position info
        
        // Calculate geometric properties
        if (!calculateShowerGeometry(shower_trunk, showercluster, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate geometry for shower " << i << std::endl;
            continue;
        }
        
        // Calculate angular properties
        if (!calculateShowerAngles(shower_trunk, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate angles for shower " << i << std::endl;
            continue;
        }
        
        // Calculate charge
        if (!calculateShowerCharge(larcv_io, event_data, reco_data, vtxIdx, i)) {
            std::cerr << "ShowerProcessor: Failed to calculate charge for shower " << i << std::endl;
            // Continue with defaults
        }
        
        // Calculate energy from shower data
        if (!calculateShowerEnergy(nuvtx, i, event_data)) {
            std::cerr << "ShowerProcessor: Failed to calculate energy for shower " << i << std::endl;
            continue;
        }
        
        // Run LArPID for particle classification if available
        int num_good_planes = 0;
        if (larpid_cnn_ != nullptr) {
            if (!runLArPID(larcv_io, shower_trunk, showercluster, i, num_good_planes, event_data, reco_data)) {
                // Fall back to default values if CNN fails
                setDefaultPIDScores(i, event_data);
            }
        } else {
            // Set default values if no CNN available
            setDefaultPIDScores(i, event_data);
        }
        event_data->showerNGoodPlanes[i] = num_good_planes;
        
        // Update energy based on PID classification
        updateEnergyBasedOnPID(nuvtx, i, event_data);
        
        // Secondary shower info
        if (i < (int)nuvtx.shower_isSecondary_v.size()) {
            event_data->showerIsSecondary[i] = nuvtx.shower_isSecondary_v.at(i);
        } else {
            event_data->showerIsSecondary[i] = isSecondaryShower(shower_trunk, showercluster);
        }
        
        // Calculate charge fractions
        calculateChargeFractions(larcv_io, i, event_data, reco_data);
        
        // Number of hits
        event_data->showerNHits[i] = (int)showercluster.size();
        
    }
    
    return true;
}

bool ShowerProcessor::calculateShowerGeometry(const larlite::track& shower_trunk,
                                             const larlite::larflowcluster& cluster,
                                             int shower_idx,
                                             EventData* event_data) {
    
    // Check containment
    bool iscontained = true;
    for (size_t ihit = 0; ihit < cluster.size(); ihit++) {
        auto const& hit = cluster[ihit];
        if (!WCFiducial::getME()->insideFV(hit[0], hit[1], hit[2])) {
            iscontained = false;
            break;
        }
    }
    event_data->showerIsContainedInFV[shower_idx] = iscontained ? 1 : 0;
    
    // Start position from shower trunk
    auto start_pos = shower_trunk.Vertex();
    event_data->showerStartPosX[shower_idx] = start_pos.X();
    event_data->showerStartPosY[shower_idx] = start_pos.Y();
    event_data->showerStartPosZ[shower_idx] = start_pos.Z();
    
    // Start direction from shower trunk
    auto start_dir = shower_trunk.VertexDirection();
    event_data->showerStartDirX[shower_idx] = start_dir.X();
    event_data->showerStartDirY[shower_idx] = start_dir.Y();
    event_data->showerStartDirZ[shower_idx] = start_dir.Z();
    
    // Distance to vertex
    event_data->showerDistToVtx[shower_idx] = calculateDistanceToVertex(shower_trunk);
    
    return true;
}

bool ShowerProcessor::calculateShowerAngles(const larlite::track& shower_trunk, int shower_idx,
                                           EventData* event_data) {
    
    // Calculate cosine of angle with beam direction (z-axis)
    event_data->showerCosTheta[shower_idx] = calculateCosTheta(shower_trunk);
    
    // Calculate cosine of angle with gravity direction (negative y-axis)
    event_data->showerCosThetaY[shower_idx] = calculateCosThetaY(shower_trunk);
    
    return true;
}

bool ShowerProcessor::calculateShowerCharge(larcv::IOManager* larcv_io,
                                           EventData* event_data,
                                           RecoData* reco_data,
                                           int vtxIdx,
                                           int shower_idx) {
    
    auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
    auto const& cluster = nuvtx.shower_v.at(shower_idx);

    auto ev_img = (larcv::EventImage2D*)larcv_io->get_data("image2d", "wire");
    auto image2Dvec = ev_img->as_vector();

    float clusterCharge = 0.0f;
    float threshold = 10.0f; // ADC threshold
    std::set< std::array<int,3> > pixels_visited; // prevent double counting pixels
    for (auto const& hit : cluster) {
        for (int p = 0; p < 3; ++p) {
            auto const& img = image2Dvec.at(p);
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
    }

    event_data->showerCharge[shower_idx] = clusterCharge;

    return true;
}

bool ShowerProcessor::calculateShowerEnergy(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx,
                                           EventData* event_data) {
    
    // Calculate different energy estimators
    float electron_energy = calculateElectronEnergy(nuvtx, shower_idx);
    float photon_energy = calculatePhotonEnergy(nuvtx, shower_idx);
    
    // Per-plane energy reconstruction
    float energy_u = calculateShowerEnergyByPlane(nuvtx, shower_idx, 0); // U plane
    float energy_v = calculateShowerEnergyByPlane(nuvtx, shower_idx, 1); // V plane  
    float energy_y = calculateShowerEnergyByPlane(nuvtx, shower_idx, 2); // Y plane
    
    // Store energies
    event_data->showerRecoE[shower_idx] = electron_energy; // Default to electron energy
    event_data->showerRecoEU[shower_idx] = energy_u;
    event_data->showerRecoEV[shower_idx] = energy_v;
    event_data->showerRecoEY[shower_idx] = energy_y;
    
    // Store shower length (approximation)
    event_data->showerLength[shower_idx] = calculateShowerLength(nuvtx, shower_idx);
    
    // Store opening angle if available
    event_data->showerOpeningAngle[shower_idx] = calculateOpeningAngle(nuvtx, shower_idx);
    
    return true;
}

bool ShowerProcessor::runLArPID(larcv::IOManager* larcv_io,
                               const larlite::track& shower_trunk,
                               const larlite::larflowcluster& cluster,
                               int shower_idx,
                               int& num_good_planes,
                               EventData* event_data,
                               RecoData* reco_data) {
    
    num_good_planes = 0;
    
    if (!larpid_cnn_) return false;

    try {
        // Get vertex information
        int vtxIdx = event_data->vtxIndex;
        auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
        
        // Get shower start point from trunk
        auto start_pos = shower_trunk.Vertex();
        
        // Run LArPID
        bool preserve_shower_pixels = true; // For showers, we want to preserve shower pixels
        bool success = false;   

        // Get spacepoints
        std::vector< std::vector<float> > hitcluster;
        hitcluster.reserve(cluster.size());
        std::cout << "shower[" << shower_idx << "] hitcluster size=" << cluster.size() << std::endl;
        for (auto const& hit : cluster) {
            std::vector<float> pt(7, 0);
            pt[0] = hit[0];
            pt[1] = hit[1];
            pt[2] = hit[2];
            pt[3] = hit.tick;
            pt[4] = hit.targetwire[0];
            pt[5] = hit.targetwire[1];
            pt[6] = hit.targetwire[2];
            hitcluster.push_back(pt);
        }

        std::vector< std::vector<larpid::data::CropPixData_t> > larpid_input
         = larpid::interface::make_prongCNN_input_sparse_images(*larcv_io, 
            hitcluster, start_pos, preserve_shower_pixels);
        
        int ngood_planes = 0;
        for (size_t p = 0; p < 3; p++) {
            std::cout << "shower[" << shower_idx << "] plane[" << p << "] npixels=" << larpid_input[p].size() << std::endl;
            if (larpid_input[p].size() >= 10)
                ngood_planes++;
        }
        std::cout << "shower[" << shower_idx << "] num good planes=" << ngood_planes << std::endl;
        num_good_planes = ngood_planes;

        larpid::data::ModelOutput output;
        
        if (ngood_planes >= 2) {
            try {
                output = larpid_cnn_->run_inference(larpid_input);
                success = true;
                std::cout << "LArPID inference successful for shower" << std::endl;
            }
            catch (std::exception& e) {
                success = false;
                std::stringstream errmsg;
                errmsg << "Error running LArPID model on shower" << std::endl;
                errmsg << e.what() << std::endl;
                throw std::runtime_error(errmsg.str());
            }
        }

        if ( is_mc_ ) {
            // get prong larpid groundtruth
            getMCProngParticles(larcv_io, larpid_input, event_data, shower_idx );
        }
        
        if (success) {
            // Assign scores
            event_data->showerElScore[shower_idx] = output.classScores[0];
            event_data->showerPhScore[shower_idx] = output.classScores[1];
            event_data->showerMuScore[shower_idx] = output.classScores[2];
            event_data->showerPiScore[shower_idx] = output.classScores[3];
            event_data->showerPrScore[shower_idx] = output.classScores[4];
            
            // Determine PID
            event_data->showerPID[shower_idx] = output.predictedPID;
            
            // Reco Quality metrics
            event_data->showerComp[shower_idx] = output.completeness;
            event_data->showerPurity[shower_idx] = output.purity;
            
            // Process type scores
            event_data->showerProcess[shower_idx] = output.predictedProcess;
            event_data->showerPrimaryScore[shower_idx] = output.processScores[0];
            event_data->showerFromChargedScore[shower_idx] = output.processScores[2];
            event_data->showerFromNeutralScore[shower_idx] = output.processScores[1];

            // Mark as classified
            event_data->showerClassified[shower_idx] = 1;
            
            return true;
        }
        else {
            event_data->showerClassified[shower_idx] = 0;
        }
    
    } catch (const std::exception& e) {
        std::cerr << "ShowerProcessor: LArPID failed with error: " 
                  << e.what() << std::endl;
    }
    
    return false;
}

void ShowerProcessor::setDefaultPIDScores(int shower_idx, EventData* event_data) {
    // Default values for unclassified showers
    event_data->showerElScore[shower_idx] = -99.0f;
    event_data->showerPhScore[shower_idx] = -99.0f;
    event_data->showerMuScore[shower_idx] = -99.0f;
    event_data->showerPiScore[shower_idx] = -99.0f;
    event_data->showerPrScore[shower_idx] = -99.0f;
    event_data->showerPID[shower_idx] = -1;
    
    // Default quality metrics
    event_data->showerComp[shower_idx] = -1.0f;
    event_data->showerPurity[shower_idx] = -1.0f;
    event_data->showerProcess[shower_idx] = -1;
    
    // Default origin scores
    event_data->showerPrimaryScore[shower_idx] = -99.0f;
    event_data->showerFromNeutralScore[shower_idx] = -99.0f;
    event_data->showerFromChargedScore[shower_idx] = -99.0f;
    
    // Mark as unclassified
    event_data->showerClassified[shower_idx] = 0;
}

void ShowerProcessor::updateEnergyBasedOnPID(const larflow::reco::NuVertexCandidate& nuvtx,
                                           int shower_idx,
                                           EventData* event_data) {
    
    int pid = event_data->showerPID[shower_idx];
    
    // For showers, energy reconstruction is mainly calorimetric
    // PID mainly affects which calibration to use
    switch (pid) {
        case PID_ELECTRON:
            // Use electron energy calculation (already stored)
            break;
        case PID_PHOTON:
            // Apply photon-specific correction
            event_data->showerRecoE[shower_idx] = calculatePhotonEnergy(nuvtx, shower_idx);
            break;
        default:
            // Keep electron energy as default for other particles
            break;
    }
}

void ShowerProcessor::calculateChargeFractions(larcv::IOManager* larcv_io,
                                             int shower_idx,
                                             EventData* event_data,
                                             RecoData* reco_data) {
    
    // Get total charge in the event (simplified)
    float total_charge = 0.0f;
    for (int i = 0; i < event_data->nShowers; i++) {
        if (event_data->showerCharge[i] > 0) {
            total_charge += event_data->showerCharge[i];
        }
    }
    
    // Calculate charge fraction
    if (total_charge > 0) {
        event_data->showerChargeFrac[shower_idx] = 
            event_data->showerCharge[shower_idx] / total_charge;
    } else {
        event_data->showerChargeFrac[shower_idx] = 0.0f;
    }
    
    // Hit fraction (placeholder)
    event_data->showerHitFrac[shower_idx] = 0.15f;
}

int ShowerProcessor::getPIDFromScores(float el_score, float ph_score, 
                                    float mu_score, float pi_score, 
                                    float pr_score) const {
    
    // Find highest score
    float max_score = el_score;
    int pid = PID_ELECTRON;
    
    if (ph_score > max_score) {
        max_score = ph_score;
        pid = PID_PHOTON;
    }
    if (mu_score > max_score) {
        max_score = mu_score;
        pid = 13; // Muon
    }
    if (pi_score > max_score) {
        max_score = pi_score;
        pid = 211; // Pion
    }
    if (pr_score > max_score) {
        max_score = pr_score;
        pid = 2212; // Proton
    }
    
    return pid;
}

float ShowerProcessor::calculateDistanceToVertex(const larlite::track& shower_trunk) const {
    auto start_pos = shower_trunk.Vertex();
    
    float dx = start_pos.X() - vertex_x_;
    float dy = start_pos.Y() - vertex_y_;
    float dz = start_pos.Z() - vertex_z_;
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

float ShowerProcessor::calculateCosTheta(const larlite::track& shower_trunk) const {
    auto direction = shower_trunk.VertexDirection();
    
    // Dot product with beam direction (0, 0, 1)
    float cos_theta = direction.Z() / direction.Mag();
    
    return cos_theta;
}

float ShowerProcessor::calculateCosThetaY(const larlite::track& shower_trunk) const {
    auto direction = shower_trunk.VertexDirection();
    
    // Dot product with gravity direction (0, -1, 0)
    float cos_theta_y = -direction.Y() / direction.Mag();
    
    return cos_theta_y;
}

std::vector<float> ShowerProcessor::getShowerDirection(const larlite::track& shower_trunk) const {
    auto direction = shower_trunk.VertexDirection();
    
    return {static_cast<float>(direction.X()), 
            static_cast<float>(direction.Y()), 
            static_cast<float>(direction.Z())};
}

float ShowerProcessor::calculateElectronEnergy(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const {
    // Get shower energy from NuVertexCandidate shower momentum vectors
    // This comes from the calorimetric reconstruction in larflow
    
    // Try to get energy from shower_plane_mom_vv if available
    if (shower_idx < (int)nuvtx.shower_plane_mom_vv.size()) {
        const auto& plane_mom_v = nuvtx.shower_plane_mom_vv[shower_idx];
        
        // Sum energy across all planes, or use the best plane (collection plane = plane 2)
        if (plane_mom_v.size() > 2) {
            // Use collection plane (Y-plane) energy if available
            return plane_mom_v[2].E(); // Energy in MeV
        } else if (plane_mom_v.size() > 0) {
            // Fall back to first available plane energy
            return plane_mom_v[0].E();
        }
    }
    
    // Fall back to calculating from cluster charge if available
    if (shower_idx < (int)nuvtx.shower_v.size()) {
        const auto& cluster = nuvtx.shower_v[shower_idx];
        float total_charge = 0.0f;
        for (const auto& hit : cluster) {
            // Simple charge summation - would need proper calibration
            total_charge += 1.0f; // placeholder charge per hit
        }
        // Convert charge to energy using approximate conversion factor
        return total_charge * 0.5f; // MeV per unit charge (rough approximation)
    }
    
    return 100.0f; // Default fallback energy in MeV
}

float ShowerProcessor::calculatePhotonEnergy(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const {
    // For photons, energy calculation is similar to electrons
    // but may have different calibration factors
    
    float electron_energy = calculateElectronEnergy(nuvtx, shower_idx);
    
    // Apply photon-specific correction factor (if any)
    // This would be calibrated based on MC studies
    float photon_correction = 1.0f; // No correction for now
    
    return electron_energy * photon_correction;
}

float ShowerProcessor::calculateShowerEnergyByPlane(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx, int plane) const {
    // Get energy for specific wire plane from NuVertexCandidate
    
    // Use shower_plane_mom_vv to get plane-specific energy
    if (shower_idx < (int)nuvtx.shower_plane_mom_vv.size()) {
        const auto& plane_mom_v = nuvtx.shower_plane_mom_vv[shower_idx];
        
        if (plane >= 0 && plane < (int)plane_mom_v.size()) {
            return plane_mom_v[plane].E(); // Energy in MeV for this plane
        }
    }
    
    // Fallback to total energy distribution if plane-specific data not available
    float total_energy = calculateElectronEnergy(nuvtx, shower_idx);
    
    // Simple approximation: distribute energy based on plane
    // In practice, this would use plane-specific calibrations
    if (plane == 2) { // Collection plane (Y) typically has most energy
        return total_energy * 0.5f;
    } else { // Induction planes (U, V)
        return total_energy * 0.25f;
    }
}

float ShowerProcessor::calculateShowerLength(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const {
    // Approximate shower length calculation
    float energy = calculateElectronEnergy(nuvtx, shower_idx);
    
    // Rough approximation: shower length scales with energy
    // For electrons in LAr: L ≈ 0.3 * E^0.5 (cm for E in MeV)
    return 0.3f * std::sqrt(energy);
}

float ShowerProcessor::calculateOpeningAngle(const larflow::reco::NuVertexCandidate& nuvtx, int shower_idx) const {
    // Placeholder: typical EM shower opening angle
    return 0.1f; // radians, approximately 6 degrees
}

bool ShowerProcessor::isSecondaryShower(const larlite::track& shower_trunk, const larlite::larflowcluster& cluster) const {
    // Simple heuristic: showers starting far from vertex are likely secondary
    float dist_to_vertex = calculateDistanceToVertex(shower_trunk);
    
    // Also consider cluster size threshold - very small clusters might be secondary
    int cluster_size = cluster.size();
    
    return (dist_to_vertex > 3.0f) || (cluster_size < 10); // 3cm distance or 10 hits threshold
}

bool ShowerProcessor::getMCProngParticles( larcv::IOManager* larcv_io,
    std::vector< std::vector<larpid::data::CropPixData_t> >& prong_vv,
    EventData* event_data,
    int shower_idx )
{

    struct ShowerInfo_t {
        int pdg;
        int nodeidx;
        float pixI;
        ShowerInfo_t()
        : pdg(-1), nodeidx(-1), pixI(0.0)
        {};
    };  
    float totalPixI = 0.0;
    std::map< int, float > particleDict;
    std::map< int, ShowerInfo_t > ShowerDict;

    int nsparse_imgs = prong_vv.size();
    //   print("[getMCProngParticle] num sparse images=",sparseimg_vv.size(),flush=True)
    //   print("  adc_v.size()=",adc_v.size(),flush=True)
    //   print("  pmcpg: ",mcpg,flush=True)
    //   print("  pmcpm: ",mcpm,flush=True)

    for ( int p=0; p<3; p++ ) { 
        auto& sparseimg = prong_vv.at(p);
        int npix = (int)sparseimg.size();
        for (int iipix=0; iipix<npix; iipix++) {
            auto& pix = sparseimg.at(iipix);
            totalPixI += pix.val;
            auto pixContents = _mcpm->getPixContent(p, pix.rawRow, pix.rawCol);   
            for ( auto& part : pixContents.particles) {
                int pdg = abs(part.pdg);
                auto it_particle = particleDict.find( pdg );
                if ( it_particle==particleDict.end() ) {
                    particleDict[ pdg ] = 0.;
                }
                particleDict[ pdg ] += pix.val;    

                auto it_Showerid = ShowerDict.find( part.tid );
                if ( it_Showerid==ShowerDict.end() ) {
                    ShowerDict[part.tid] = ShowerInfo_t();
                    ShowerDict[part.tid].pdg = part.pdg;
                    ShowerDict[part.tid].nodeidx = part.nodeidx;
                }
                ShowerDict[part.tid].pixI += pixContents.pixI;
            }
        }
    }

    int maxPartPDG = 0; 
    int maxPartNID = -1;
    int maxPartTID = -1;
    int maxPartMID = -1;
    float maxPartI    = 0.;
    float maxPartComp = 0.;
    float maxPartE    = -1.;
    std::vector<int> pdglist;
    std::vector<int> puritylist;

    for (auto it_part=particleDict.begin(); it_part!=particleDict.end(); it_part++ ) {
        pdglist.push_back( it_part->first );
        float purity = 0.;
        if ( totalPixI>0.0 ) {
            purity = (it_part->second)/totalPixI;
        }
        puritylist.push_back(purity);
    }

    for ( auto it_Shower=ShowerDict.begin(); it_Shower!=ShowerDict.end(); it_Shower++ ) {
        std::cout << "Shower[" << shower_idx << "] pdg=" << it_Shower->second.pdg << " totalpix=" << it_Shower->second.pixI << std::endl; 
        if ( it_Shower->second.pixI > maxPartI ) {
            maxPartI   = it_Shower->second.pixI;
            maxPartPDG = it_Shower->second.pdg;
            maxPartNID = it_Shower->second.nodeidx;
            maxPartTID = it_Shower->first;
        }
    }

    float totNodePixI = 0.;
    auto ev_adc = (larcv::EventImage2D*)larcv_io->get_data(larcv::kProductImage2D,"wire");
    auto& adc_v = ev_adc->as_vector();
    if ( maxPartI>0. ) {
        auto& maxPartNode = _mcpg->node_v.at(maxPartNID);

        maxPartMID = maxPartNode.mtid;
        maxPartE = maxPartNode.E_MeV;
        if ( maxPartNode.tid != maxPartTID ) {
            throw std::runtime_error( "ERROR: mismatch between node Shower id from mcpm and mcpg in getMCProngParticle" );
        }
        for (int p=0; p<3; p++) {
            auto& pixels = maxPartNode.pix_vv.at(p);
            for (int iP=0; iP<(int)pixels.size()/2; iP++ ) {
                int row = ( pixels[2*iP]-2400 )/6;
                int col = pixels[2*iP+1];
                totNodePixI += adc_v[p].pixel(row,col);
            }
        }
        if ( totNodePixI>0) {
            maxPartComp = maxPartI/totNodePixI;
        }
        if ( maxPartComp>1.0) {
            std::cout << "ERROR: prong completeness calculated to be >1" << std::endl;
        }
    }
    float maxPartPurity = (totalPixI) ? maxPartI/totalPixI : 0.0;

    event_data->showerTruePID[shower_idx]    = maxPartPDG;
    event_data->showerTrueTID[shower_idx]    = maxPartTID;
    event_data->showerTrueMID[shower_idx]    = maxPartMID;
    event_data->showerTrueE[shower_idx]      = maxPartE;
    event_data->showerTruePurity[shower_idx] = maxPartPurity;
    event_data->showerTrueComp[shower_idx]   = maxPartComp;
    
    event_data->showerTrueElPurity[shower_idx] = 0.;
    event_data->showerTruePhPurity[shower_idx] = 0.;
    event_data->showerTrueMuPurity[shower_idx] = 0.;
    event_data->showerTruePiPurity[shower_idx] = 0.;
    event_data->showerTruePrPurity[shower_idx] = 0.;

    int ii=0;
    for ( auto& pdg : pdglist ) {
        switch (pdg) {
        case 11:
            event_data->showerTrueElPurity[shower_idx] = puritylist.at(ii);
            break;
        case 22:
            event_data->showerTruePhPurity[shower_idx] = puritylist.at(ii);
            break;
        case 13:
            event_data->showerTrueMuPurity[shower_idx] = puritylist.at(ii);
            break;
        case 211:
            event_data->showerTruePiPurity[shower_idx] = puritylist.at(ii);
            break;
        case 2212:
            event_data->showerTruePrPurity[shower_idx] = puritylist.at(ii);
            break;         
        }
    }

    return true;

}

} // namespace gen2ntuple