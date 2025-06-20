#include "VertexSummaryProcessor.h"

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
//       #project event 3D points to first PC axis with 0 at vertex projection
//       vtxProj = [0.,0.,0.]
//       ldist = 0.
//       for c in range(3):
//         ldist += (vertex.pos[c] - eventCluster.pca_center[c])*eventPCAxis0[c]
//       for c in range(3):
//         vtxProj[c] = eventCluster.pca_center[c] + ldist*eventPCAxis0[c]

//       projDists = []
//       avgPosPDistTimes = 0
//       avgNegPDistTimes = 0
//       nPosPDists = 0
//       nNegPDists = 0
//       for hit in eventLarflowCluster:
//         projPt = [0.,0.,0.]
//         ldist = 0.
//         for c in range(3):
//           ldist += (hit[c] - eventCluster.pca_center[c])*eventPCAxis0[c]
//         if ldist > 0.:
//           avgPosPDistTimes += hit.tick
//           nPosPDists += 1
//         else:
//           avgNegPDistTimes += hit.tick
//           nNegPDists += 1
//         for c in range(3):
//           projPt[c] = eventCluster.pca_center[c] + ldist*eventPCAxis0[c]
//         projDist = 0.
//         for c in range(3):
//           projDist += (projPt[c] - vtxProj[c])**2
//         projDists.append(sqrt(projDist))
//       projDists.sort()

//       #check if PCA axis is pointing in direction of increaasing hit times
//       avgPosPDistTimes /= (1.0*nPosPDists)
//       avgNegPDistTimes /= (1.0*nNegPDists)
//       if avgPosPDistTimes > 0.:
//         eventPCAxis0TSlope[0] = 1
//       else:
//         eventPCAxis0TSlope[0] = -1

//       #calculate maximum point gap along PCA projection and maximum charge in between gaps
//       maxGapFull = -1.
//       maxGap90 = -1.
//       maxGap80 = -1.
//       maxGap70 = -1.
//       maxGap60 = -1.
//       maxCntD02 = -1.
//       maxCntD04 = -1.
//       maxCntD06 = -1.
//       maxCntD08 = -1.
//       maxCntD10 = -1.
//       currentCntD02 = 0.
//       currentCntD04 = 0.
//       currentCntD06 = 0.
//       currentCntD08 = 0.
//       currentCntD10 = 0.

//       for iEP in range(1, len(projDists)):

//         gap = projDists[iEP] - projDists[iEP-1]

//         if gap > maxGapFull:
//           maxGapFull = gap
//         if iEP < int(0.9*len(projDists)) and gap > maxGap90:
//           maxGap90 = gap
//         if iEP < int(0.8*len(projDists)) and gap > maxGap80:
//           maxGap80 = gap
//         if iEP < int(0.7*len(projDists)) and gap > maxGap70:
//           maxGap70 = gap
//         if iEP < int(0.6*len(projDists)) and gap > maxGap60:
//           maxGap60 = gap

//         if gap > 2. or iEP == (len(projDists) - 1):
//           if currentCntD02 > maxCntD02:
//             maxCntD02 = currentCntD02
//           currentCntD02 = 0.
//         else:
//           currentCntD02 += gap
//         if gap > 4. or iEP == (len(projDists) - 1):
//           if currentCntD04 > maxCntD04:
//             maxCntD04 = currentCntD04
//           currentCntD04 = 0.
//         else:
//           currentCntD04 += gap
//         if gap > 6. or iEP == (len(projDists) - 1):
//           if currentCntD06 > maxCntD06:
//             maxCntD06 = currentCntD06
//           currentCntD06 = 0.
//         else:
//           currentCntD06 += gap
//         if gap > 8. or iEP == (len(projDists) - 1):
//           if currentCntD08 > maxCntD08:
//             maxCntD08 = currentCntD08
//           currentCntD08 = 0.
//         else:
//           currentCntD08 += gap
//         if gap > 10. or iEP == (len(projDists) - 1):
//           if currentCntD10 > maxCntD10:
//             maxCntD10 = currentCntD10
//           currentCntD10 = 0.
//         else:
//           currentCntD10 += gap

//       eventPCProjMaxGap[0] = maxGapFull
//       eventPCProjMaxGap[1] = maxGap90
//       eventPCProjMaxGap[2] = maxGap80
//       eventPCProjMaxGap[3] = maxGap70
//       eventPCProjMaxGap[4] = maxGap60
//       eventPCProjMaxDist[0] = maxCntD02
//       eventPCProjMaxDist[1] = maxCntD04
//       eventPCProjMaxDist[2] = maxCntD06
//       eventPCProjMaxDist[3] = maxCntD08
//       eventPCProjMaxDist[4] = maxCntD10

    return true;

}


}