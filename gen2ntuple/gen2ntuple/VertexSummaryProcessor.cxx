#include "VertexSummaryProcessor.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "larflow/RecoUtils/cluster_t.h"
#include "larflow/RecoUtils/cluster_functions.h"

namespace gen2ntuple {

/**
 * @brief Calculate hit and charge fraction per shower and track 
 * 
 * Assumes we've already processed Track and Showers
 *
 */
bool VertexSummaryProcessor::calculateHitAndChargeFractions( EventData* event_data ) 
{

    // Get total charge in the event
    float total_charge = 0.0f;
    float total_hit = 0;

    for (int i = 0; i < event_data->nTracks; i++) {
        if (event_data->trackCharge[i] > 0) {
            total_charge += event_data->trackCharge[i];
        }
        if (event_data->trackNHits[i]>0) {
            total_hit += (float)event_data->trackNHits[i];
        }
    }

    for (int i=0; i < event_data->nShowers; i++) {
        if (event_data->showerCharge[i]>0) {
            total_charge += event_data->showerCharge[i];
        }
        if (event_data->showerNHits[i]>0) {
            total_hit += (float)event_data->showerNHits[i];
        }
    }
    
    // Calculate fractions: tracks
    for (int i=0;  i < event_data->nTracks; i++ ) {
        if (total_charge > 0) {
            event_data->trackChargeFrac[i] = 
                event_data->trackCharge[i] / total_charge;
        }
        else{
            event_data->trackChargeFrac[i] = 0.0f;
        }

        if ( total_hit>0 ) {
            event_data->trackHitFrac[i] =
                (float)event_data->trackNHits[i] / total_hit;
        } else {
            event_data->trackHitFrac[i] = 0.0f;
        }
    }   

    // Calculate fractions: showers
    for (int i=0;  i < event_data->nShowers; i++ ) {
        if (total_charge > 0) {
            event_data->showerChargeFrac[i] = 
                event_data->showerCharge[i] / total_charge;
        }
        else{
            event_data->showerChargeFrac[i] = 0.0f;
        }

        if ( total_hit>0 ) {
            event_data->showerHitFrac[i] =
                (float)event_data->showerNHits[i] / total_hit;
        } else {
            event_data->showerHitFrac[i] = 0.0f;
        }
    }

    return true;

}

bool VertexSummaryProcessor::calculateEventPCA( EventData* event_data,
        RecoData* reco_data )
{

    if ( event_data->vtxIndex < 0) {
        // no vertex. return.
        return true;
    }

    auto const& nuvtx = reco_data->nuvtx_v->at( event_data->vtxIndex );

    int total_hits = 0;
    for (int i=0; i<(int)nuvtx.track_hitcluster_v.size(); i++) {
        total_hits += nuvtx.track_hitcluster_v[i].size();
    }
    for (int i=0; i<(int)nuvtx.shower_v.size(); i++) {
        total_hits += nuvtx.shower_v[i].size();
    }

    larflow::recoutils::cluster_t eventCluster;
    eventCluster.points_v.reserve( total_hits );

    for (int i=0; i<(int)nuvtx.track_hitcluster_v.size(); i++) {
        for ( auto& hit : nuvtx.track_hitcluster_v.at(i) ) {
            std::vector<float> pt = { hit[0], hit[1], hit[2] };
            eventCluster.points_v.push_back(pt);
        }
    }
    for (int i=0; i<(int)nuvtx.shower_v.size(); i++) {
        for ( auto& hit : nuvtx.shower_v.at(i) ) {
            std::vector<float> pt = { hit[0], hit[1], hit[2] };
            eventCluster.points_v.push_back(pt);
        }        
    }

    try {
        larflow::recoutils::cluster_pca( eventCluster );
    }
    catch ( std::exception& e ) {
        std::cerr << "Could not performce PCA on event cluster" << std::endl;
        return true;
    }

    for (int iPCA=0; iPCA<3; iPCA++) {
      event_data->eventPCAxis1[iPCA]     = eventCluster.pca_axis_v[1][iPCA];
      event_data->eventPCAxis0[iPCA]     = eventCluster.pca_axis_v[0][iPCA];
      event_data->eventPCAxis2[iPCA]     = eventCluster.pca_axis_v[2][iPCA];
      event_data->eventPCEigenVals[iPCA] = eventCluster.pca_eigenvalues[iPCA];
    }
    
    // Project event 3D points to first PC axis with 0 at vertex projection
    std::vector<float> vtxProj(3, 0.0f);
    float ldist = 0.0f;
    // Calculate projection distance from vertex to PCA center along first PC axis
    ldist += (event_data->vtxX - eventCluster.pca_center[0]) * event_data->eventPCAxis0[0];
    ldist += (event_data->vtxY - eventCluster.pca_center[1]) * event_data->eventPCAxis0[1];
    ldist += (event_data->vtxZ - eventCluster.pca_center[2]) * event_data->eventPCAxis0[2];
    
    // Calculate vertex projection point on PC axis
    for (int c = 0; c < 3; c++) {
        vtxProj[c] = eventCluster.pca_center[c] + ldist * event_data->eventPCAxis0[c];
    }
    // Collect all projection distances
    std::vector<float> projDists;
    projDists.reserve(total_hits);
    float avePosPDistTimes = 0.f;
    float aveNegPDistTimes = 0.f;
    int   nPosPDists = 0;
    int   nNegPDists = 0;
    
    // Calculate projection distances for all hits
    for (int i = 0; i < (int)nuvtx.track_hitcluster_v.size(); i++) {
        for (auto& hit : nuvtx.track_hitcluster_v.at(i)) {
            std::vector<float> projPt(3, 0.0f);
            float ldist = 0.0f;
            for (int c = 0; c < 3; c++) {
                ldist += (hit[c] - eventCluster.pca_center[c]) * event_data->eventPCAxis0[c];
            }
            if ( ldist>0 ) {
                avePosPDistTimes += hit.tick;
                nPosPDists++;
            }
            else {
                aveNegPDistTimes += hit.tick;
                nNegPDists += 1;
            }

            for (int c = 0; c < 3; c++) {
                projPt[c] = eventCluster.pca_center[c] + ldist * event_data->eventPCAxis0[c];
            }
            float projDist = 0.0f;
            for (int c = 0; c < 3; c++) {
                projDist += (projPt[c] - vtxProj[c]) * (projPt[c] - vtxProj[c]);
            }
            projDists.push_back(std::sqrt(projDist));
        }
    }
    
    for (int i = 0; i < (int)nuvtx.shower_v.size(); i++) {
        for (auto& hit : nuvtx.shower_v.at(i)) {
            std::vector<float> projPt(3, 0.0f);
            float ldist = 0.0f;
            for (int c = 0; c < 3; c++) {
                ldist += (hit[c] - eventCluster.pca_center[c]) * event_data->eventPCAxis0[c];
            }

            if ( ldist>0 ) {
                avePosPDistTimes += hit.tick;
                nPosPDists++;
            }
            else {
                aveNegPDistTimes += hit.tick;
                nNegPDists += 1;
            }

            for (int c = 0; c < 3; c++) {
                projPt[c] = eventCluster.pca_center[c] + ldist * event_data->eventPCAxis0[c];
            }
            float projDist = 0.0f;
            for (int c = 0; c < 3; c++) {
                projDist += (projPt[c] - vtxProj[c]) * (projPt[c] - vtxProj[c]);
            }
            projDists.push_back(std::sqrt(projDist));
        }
    }
    
    // Sort projection distances
    std::sort(projDists.begin(), projDists.end());

    // #check if PCA axis is pointing in direction of increaasing hit times
    avePosPDistTimes /= (1.0*nPosPDists);
    aveNegPDistTimes /= (1.0*nNegPDists);
    if ( avePosPDistTimes > 0.)
      event_data->eventPCAxis0TSlope = 1;
    else
      event_data->eventPCAxis0TSlope = -1;

    // Calculate maximum point gap along PCA projection
    float maxGapFull = -1.0f;
    float maxCntD02 = -1.0f;
    float currentCntD02 = 0.0f;
    
    for (size_t iEP = 1; iEP < projDists.size(); iEP++) {
        float gap = projDists[iEP] - projDists[iEP-1];
        
        // Track maximum gap
        if (gap > maxGapFull) {
            maxGapFull = gap;
        }
        
        // Track maximum continuous distance with gaps <= 2.0
        if (gap > 2.0f || iEP == (projDists.size() - 1)) {
            if (currentCntD02 > maxCntD02) {
                maxCntD02 = currentCntD02;
            }
            currentCntD02 = 0.0f;
        } else {
            currentCntD02 += gap;
        }
    }
    
    // Store the results
    event_data->eventPCProjMaxGap = maxGapFull;
    event_data->eventPCProjMaxDist = maxCntD02;

    return true;

}


}