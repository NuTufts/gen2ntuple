#include "TrackProcessor.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <array>
#include <set>

// LArLite includes
#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/track.h"
#include "larlite/DataFormat/vertex.h"
#include "larlite/DataFormat/mctrack.h"
#include "larlite/DataFormat/mcpart.h"
#include "larlite/DataFormat/larflowcluster.h"

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
#include "pionRange2T.h"

namespace gen2ntuple {

TrackProcessor::TrackProcessor() 
    : is_mc_(false), vertex_x_(0), vertex_y_(0), vertex_z_(0), larpid_cnn_(nullptr) {
}

bool TrackProcessor::processEvent( larlite::storage_manager* larlite_io,
                                   larcv::IOManager* larcv_io,
                                   EventData* event_data,
                                   RecoData*  reco_data ) {
    
    if (!event_data) {
        std::cerr << "TrackProcessor: EventData pointer is null" << std::endl;
        return false;
    }
    
    // Set vertex position from event data
    setVertexPosition(event_data->vtxX, event_data->vtxY, event_data->vtxZ);
    
    // Extract track information from LArLite
    if (!extractTrackInfo(larlite_io, larcv_io, event_data, reco_data)) {
        // No tracks found - set defaults
        event_data->nTracks = 0;
        return true;
    }
    
    return true;
}

bool TrackProcessor::extractTrackInfo( larlite::storage_manager* larlite_io,
                                       larcv::IOManager* larcv_io,
                                       EventData* event_data,
                                       RecoData*  reco_data ) 
{
    
    // Get track collection from LArLite
    int vtxIdx = event_data->vtxIndex;

    if ( reco_data->nuvtx_v->size()==0 ) {
        // no reco vertices in this event. just return.
        std::cout << "TrackProcessor::extractTrackInfo - no vertex in event" << std::endl;
        return true;
    }

    if ( vtxIdx<0 || vtxIdx>=(int)reco_data->nuvtx_v->size() ) {
        throw std::runtime_error("TrackProcess::extractTrackInfo -- vtxIndex unset -- needs to be selected in VertexSelector.");
    }
    
    auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
    
    int n_tracks = std::min((int)nuvtx.track_v.size(), EventData::MAX_TRACKS);
    event_data->nTracks = n_tracks;

    std::cout << "TrackProcessor::extractTrackInfo - number of Tracks = " << n_tracks << std::endl;
    
    // Process each track
    for (int i = 0; i < n_tracks; i++) {
        const auto& track = nuvtx.track_v.at(i);
        const auto& trackcluster = nuvtx.track_hitcluster_v.at(i);
        
        // Calculate geometric properties
        if (!calculateTrackGeometry(track, trackcluster, i, event_data)) {
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

        if (!calculateTrackCharge(larcv_io, event_data, reco_data, vtxIdx, i)) {
            std::cerr << "TrackProcessor: Failed to calculate track charge" << std::endl;
            continue;
        }
        
        // Run ProngCNN for particle classification if available
        int num_good_planes = 0;
        if (larpid_cnn_ != nullptr) {
            if (!runProngCNN(larcv_io, track, trackcluster, i, num_good_planes, event_data, reco_data)) {
                // Fall back to default values if CNN fails
                setDefaultPIDScores(i, event_data);
            }
        } else {
            // Set default values if no CNN available
            setDefaultPIDScores(i, event_data);
        }
        event_data->trackNGoodPlanes[i] = num_good_planes; 
        
        // Update energy based on PID classification
        updateEnergyBasedOnPID(vtxIdx, i, event_data, reco_data);
        
        // Secondary track info
        if (i < (int)nuvtx.track_isSecondary_v.size()) {
            event_data->trackIsSecondary[i] = nuvtx.track_isSecondary_v.at(i);
        } else {
            event_data->trackIsSecondary[i] = -1;
        }
               
    }//end of track loop

    
    return true;
}

bool TrackProcessor::calculateTrackGeometry(const larlite::track& track,
                                           const larlite::larflowcluster& cluster,
                                           int track_idx,
                                           EventData* event_data) 
{
    if (track.NumberTrajectoryPoints() < 2) {
        // Not enough points to calculate geometry
        return false;
    }

    // Check containment
    bool iscontained = true;
    for (size_t ihit = 0; ihit < cluster.size(); ihit++) {
        auto const& hit = cluster[ihit];
        if (!WCFiducial::getME()->insideFV(hit[0], hit[1], hit[2])) {
            iscontained = false;
            break;
        }
    }
    event_data->trackIsContainedInFV[track_idx] = iscontained ? 1 : 0;

    // Number of hits
    event_data->trackNHits[track_idx] = (int)cluster.size();

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
    
    // Start direction - find direction from first 5cm of track
    auto trackDirPt1 = track.Vertex();
    TVector3 trackDirPt2 = track.Vertex();
    for (int ipt = 1; ipt < (int)track.NumberTrajectoryPoints(); ipt++) {
        trackDirPt2 = track.LocationAtPoint(ipt);
        double dist = (trackDirPt2 - trackDirPt1).Mag();
        if (dist > 5.0) {
            break;
        }
    }
    auto trackDir = trackDirPt2 - trackDirPt1;
    if (trackDir.Mag() > 0.0) {
        trackDir *= (1.0 / trackDir.Mag());
    }
    event_data->trackStartDirX[track_idx] = trackDir[0];
    event_data->trackStartDirY[track_idx] = trackDir[1];
    event_data->trackStartDirZ[track_idx] = trackDir[2];
    
    // Distance to vertex
    event_data->trackDistToVtx[track_idx] = calculateDistanceToVertex(track);
    
    return true;
}

bool TrackProcessor::calculateTrackAngles(const larlite::track& track, 
                                         int track_idx,
                                         EventData* event_data) {
    
    auto direction = getTrackDirection(track);

    // Calculate cosine of angle with beam direction (z-axis)
    event_data->trackCosTheta[track_idx] = direction[2];
    
    // Calculate cosine of angle with gravity direction (negative y-axis)
    event_data->trackCosThetaY[track_idx] = -direction[1];

    return true;
}

bool TrackProcessor::calculateTrackEnergy(const larlite::track& track, int track_idx,
                                         EventData* event_data) {
    
    float track_length = calculateTrackLength(track);
    event_data->trackLength[track_idx] = track_length;
    
    // Calculate range-based energies for different particle hypotheses
    event_data->trackMuonE[track_idx] = calculateMuonEnergy(track);
    event_data->trackProtonE[track_idx] = calculateProtonEnergy(track);
    
    // Initial energy assignment (will be updated after PID)
    event_data->trackRecoE[track_idx] = event_data->trackMuonE[track_idx];
    
    return true;
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


std::vector<float> TrackProcessor::getTrackDirection(const larlite::track& track) const {
    // we follow mmr ntuple maker and use point 5 away from start

    auto trackDirPt1 = track.Vertex();
    TVector3 trackDirPt2 = track.Vertex();
    for (int ipt = 1; ipt < (int)track.NumberTrajectoryPoints(); ipt++) {
        trackDirPt2 = track.LocationAtPoint(ipt);
        double dist = (trackDirPt2 - trackDirPt1).Mag();
        if (dist > 5.0) {
            break;
        }
    }
    auto trackDir = trackDirPt2 - trackDirPt1;
    if (trackDir.Mag() > 0.0) {
        trackDir *= (1.0 / trackDir.Mag());
    }
   
    return {static_cast<float>(trackDir.X()), 
            static_cast<float>(trackDir.Y()), 
            static_cast<float>(trackDir.Z())};
}

float TrackProcessor::calculateMuonEnergy(const larlite::track& track) const {
    float length = calculateTrackLength(track);
    
    // Muon range-energy relationship from MicroBooNE
    // Based on CSDA range tables for LAr
    if (length < 0.0f) return 0.0f;
    
    // Polynomial fit for muon kinetic energy vs range
    // E = a*R + b*R^2 + c*R^3 (simplified)
    float a = 2.104f;  // MeV/cm
    float b = 0.000f;  // MeV/cm^2
    float c = 0.0022f; // MeV/cm^3
    
    float energy = a * length + b * length * length + c * length * length * length;
    return energy;
}

float TrackProcessor::calculateProtonEnergy(const larlite::track& track) const {
    float length = calculateTrackLength(track);
    
    // Proton range-energy relationship
    if (length < 0.0f) return 0.0f;
    
    // Simplified proton range-energy for LAr
    // More complex in reality due to Bragg peak
    float energy = 5.77f * std::pow(length, 0.78f);
    return energy;
}

float TrackProcessor::calculatePionRangeEnergy(float track_length) const {
    // Pion range-energy calculation
    // Based on pionRange2T function from Python
    
    if (track_length < 1.0f) return 0.0f;
    
    // Pion kinetic energy from range (simplified)
    float a = 1.953f;  // MeV/cm
    float b = 0.001f;  // MeV/cm^2
    
    float energy = a * track_length + b * track_length * track_length;
    return energy;
}

bool TrackProcessor::calculateTrackCharge(larcv::IOManager* larcv_io,
                                         EventData* event_data,
                                         RecoData* reco_data,
                                         int vtxIdx,
                                         int track_idx) 
{
    auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
    auto const& cluster = nuvtx.track_hitcluster_v.at(track_idx);

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

    event_data->trackCharge[track_idx] = clusterCharge;

    return true;
}

bool TrackProcessor::runProngCNN(larcv::IOManager* larcv_io,
                                const larlite::track& track,
                                const larlite::larflowcluster& cluster,
                                int track_idx,
                                int& num_good_planes,
                                EventData* event_data,
                                RecoData* reco_data) {
    
    num_good_planes = 0;

    if (!larpid_cnn_) return false;

    try {
        // Get vertex information
        int vtxIdx = event_data->vtxIndex;
        auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
        
        // Create prong endpoint for CNN
        std::vector<float> vtx_pos = {vertex_x_, vertex_y_, vertex_z_};
        
        // Get track start point
        auto start_pos = track.Vertex();
        std::vector<float> track_start = {
            (float)start_pos.X(), 
            (float)start_pos.Y(), 
            (float)start_pos.Z()
        };

        // Get track end point
        auto end_pos = track.LocationAtPoint( track.NumberTrajectoryPoints()-1 );
        
        // Run ProngCNN
        bool preserve_shower_pixels = false;
        bool success = false;   

        // get spacepoints
        std::vector< std::vector<float> > hitcluster;
        hitcluster.reserve( cluster.size() );
        std::cout << "track[" << track_idx << "] hitcluster size=" << cluster.size() << std::endl;
        for (auto const& hit : cluster ) {
            std::vector<float> pt(7,0);
            pt[0] = hit[0];
            pt[1] = hit[1];
            pt[2] = hit[2];
            pt[3] = hit.tick;
            pt[4] = hit.targetwire[0];
            pt[5] = hit.targetwire[1];
            pt[6] = hit.targetwire[2];
            // std::stringstream sspt;
            // sspt << "( ";
            // for ( auto const& val : pt )
            //   sspt << val << " ";
            // sspt << ")";
            // std::cout << sspt.str() << std::endl;
            hitcluster.push_back( pt );
        }

        std::vector< std::vector<larpid::data::CropPixData_t> > larpid_input
         = larpid::interface::make_prongCNN_input_sparse_images( *larcv_io, 
            hitcluster, end_pos, preserve_shower_pixels );
        
        int ngood_planes = 0;
        for (size_t p=0; p<3; p++ ) {
            std::cout << "track["  << track_idx << "] plane[" <<  p << "] npixels=" << larpid_input[p].size() << std::endl;
            if ( larpid_input[p].size()>=10 )
                ngood_planes++;
        }
        std::cout << "track[" << track_idx << "] num good planes=" << ngood_planes << std::endl;
        num_good_planes = ngood_planes;

        larpid::data::ModelOutput output;
        
        if ( ngood_planes>=2 ) {
            try {
                output = larpid_cnn_->run_inference( larpid_input );
                success = true;
                std::cout << "LArPID inference successful" << std::endl;
            }
            catch ( std::exception& e ) {
                success = false;
                std::stringstream errmsg;
                errmsg << "Error running LArPID model" << std::endl;
                errmsg << e.what() << std::endl;
                throw std::runtime_error(errmsg.str());
            }
        }
        
        if ( is_mc_ ) {
            // get prong larpid groundtruth
            getMCProngParticles(larcv_io, larpid_input, event_data, track_idx );
        }

        if (success) {
            // Assign scores
            event_data->trackElScore[track_idx] = output.classScores[0];
            event_data->trackPhScore[track_idx] = output.classScores[1];
            event_data->trackMuScore[track_idx] = output.classScores[2];
            event_data->trackPiScore[track_idx] = output.classScores[3];
            event_data->trackPrScore[track_idx] = output.classScores[4];
            
            // Determine PID
            event_data->trackPID[track_idx] = output.predictedPID;
            
            // Reco Quality metrics
            event_data->trackComp[track_idx]    = output.completeness;
            event_data->trackPurity[track_idx]  = output.purity;
            
            
            // Process type scores (if available)
            event_data->trackProcess[track_idx] = output.predictedProcess;
            event_data->trackPrimaryScore[track_idx]     = output.processScores[0];
            event_data->trackFromChargedScore[track_idx] = output.processScores[2];
            event_data->trackFromNeutralScore[track_idx] = output.processScores[1];

            // Mark as classified
            event_data->trackClassified[track_idx] = 1;
            
            return true;
        }
        else {
            event_data->trackClassified[track_idx] = 0;
        }
    
    } catch (const std::exception& e) {
        std::cerr << "TrackProcessor: ProngCNN failed with error: " 
                  << e.what() << std::endl;
    }
    
    return false;
}

void TrackProcessor::setDefaultPIDScores(int track_idx, EventData* event_data) {
    // Default to muon-like
    event_data->trackElScore[track_idx] = -99.0f;
    event_data->trackPhScore[track_idx] = -99.0f;
    event_data->trackMuScore[track_idx] = -99.0f;
    event_data->trackPiScore[track_idx] = -99.0f;
    event_data->trackPrScore[track_idx] = -99.0f;
    event_data->trackPID[track_idx]     = -1;
    
    // Default quality metrics
    event_data->trackComp[track_idx]    = -1.0f;
    event_data->trackPurity[track_idx]  = -1.0f;
    event_data->trackProcess[track_idx] = -1;
    
    // Default origin scores
    event_data->trackPrimaryScore[track_idx]     = -99.0f;
    event_data->trackFromNeutralScore[track_idx] = -99.0f;
    event_data->trackFromChargedScore[track_idx] = -99.0f;
    
    // Mark as unclassified
    event_data->trackClassified[track_idx] = 0;
}

bool TrackProcessor::updateEnergyBasedOnPID(int vtxIdx, int track_idx,
                                           EventData* event_data,
                                           RecoData* reco_data ) {
    
    auto const& nuvtx = reco_data->nuvtx_v->at(vtxIdx);
    auto const& track = nuvtx.track_v.at(track_idx);

    if ( event_data->trackPID[track_idx]<=0 ) {
        event_data->trackRecoE[track_idx] = -1.0;
        return true;
    }

    bool foundEnergy = true;
    if ( event_data->trackMuScore[track_idx] >= event_data->trackPiScore[track_idx] 
         && event_data->trackMuScore[track_idx] >= event_data->trackPrScore[track_idx] ) {
    
        event_data->trackRecoE[track_idx] = nuvtx.track_kemu_v[track_idx];
    }
    else if (event_data->trackPrScore[track_idx] >= event_data->trackMuScore[track_idx]
         && event_data->trackPrScore[track_idx] >= event_data->trackPiScore[track_idx] ) {
      
        event_data->trackRecoE[track_idx] = nuvtx.track_keproton_v[track_idx];

    }
    else {
        event_data->trackRecoE[track_idx] 
            = pionRange2T::get()->Eval(track.Length());
    }  

    return true;

}

int TrackProcessor::getPIDFromScores(float el_score, float ph_score, 
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
        pid = PID_MUON;
    }
    if (pi_score > max_score) {
        max_score = pi_score;
        pid = PID_PION;
    }
    if (pr_score > max_score) {
        max_score = pr_score;
        pid = PID_PROTON;
    }
    
    return pid;
}

bool TrackProcessor::getMCProngParticles( larcv::IOManager* larcv_io,
    std::vector< std::vector<larpid::data::CropPixData_t> >& prong_vv,
    EventData* event_data,
    int track_idx )
{

    struct TrackInfo_t {
        int pdg;
        int nodeidx;
        float pixI;
        TrackInfo_t()
        : pdg(-1), nodeidx(-1), pixI(0.0)
        {};
    };  
    float totalPixI = 0.0;
    std::map< int, float > particleDict;
    std::map< int, TrackInfo_t > trackDict;

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

                auto it_trackid = trackDict.find( part.tid );
                if ( it_trackid==trackDict.end() ) {
                    trackDict[part.tid] = TrackInfo_t();
                    trackDict[part.tid].pdg = part.pdg;
                    trackDict[part.tid].nodeidx = part.nodeidx;
                }
                trackDict[part.tid].pixI += pixContents.pixI;
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

    for ( auto it_track=trackDict.begin(); it_track!=trackDict.end(); it_track++ ) {
        std::cout << "track[" << track_idx << "] pdg=" << it_track->second.pdg << " totalpix=" << it_track->second.pixI << std::endl; 
        if ( it_track->second.pixI > maxPartI ) {
            maxPartI   = it_track->second.pixI;
            maxPartPDG = it_track->second.pdg;
            maxPartNID = it_track->second.nodeidx;
            maxPartTID = it_track->first;
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
            throw std::runtime_error( "ERROR: mismatch between node track id from mcpm and mcpg in getMCProngParticle" );
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

    event_data->trackTruePID[track_idx]    = maxPartPDG;
    event_data->trackTrueTID[track_idx]    = maxPartTID;
    event_data->trackTrueMID[track_idx]    = maxPartMID;
    event_data->trackTrueE[track_idx]      = maxPartE;
    event_data->trackTruePurity[track_idx] = maxPartPurity;
    event_data->trackTrueComp[track_idx]   = maxPartComp;
    
    event_data->trackTrueElPurity[track_idx] = 0.;
    event_data->trackTruePhPurity[track_idx] = 0.;
    event_data->trackTrueMuPurity[track_idx] = 0.;
    event_data->trackTruePiPurity[track_idx] = 0.;
    event_data->trackTruePrPurity[track_idx] = 0.;

    int ii=0;
    for ( auto& pdg : pdglist ) {
        switch (pdg) {
        case 11:
            event_data->trackTrueElPurity[track_idx] = puritylist.at(ii);
            break;
        case 22:
            event_data->trackTruePhPurity[track_idx] = puritylist.at(ii);
            break;
        case 13:
            event_data->trackTrueMuPurity[track_idx] = puritylist.at(ii);
            break;
        case 211:
            event_data->trackTruePiPurity[track_idx] = puritylist.at(ii);
            break;
        case 2212:
            event_data->trackTruePrPurity[track_idx] = puritylist.at(ii);
            break;         
        }
    }

    return true;

}


} // namespace gen2ntuple
