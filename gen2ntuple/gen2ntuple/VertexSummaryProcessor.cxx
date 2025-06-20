#include "VertexSummaryProcessor.h"

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


}