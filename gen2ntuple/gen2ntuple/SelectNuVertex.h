#ifndef __GEN2NTUPLE_SELECT_NU_VERTEX_H__
#define __GEN2NTUPLE_SELECT_NU_VERTEX_H__

#include "VertexSelectionBase.h"

namespace gen2ntuple {

    /**
     * @brief Select vertex based on highest neutrino keypoint score
     * 
     * This is the default/baseline selection method that picks the vertex
     * with the highest neutrino keypoint score (keypoint_type == 0).
     * Corresponds to the Python highest_nu_keypoint_score() function.
     */
    class SelectNuVertex : public VertexSelectionBase {
    public:
        SelectNuVertex() {}
        virtual ~SelectNuVertex() {}
        
        /**
         * @brief Select vertex with highest neutrino keypoint score
         * 
         * @param larcv_io LArCV IO manager
         * @param larlite_io LArLite storage manager
         * @param nuvtx_v Vector of neutrino vertex candidates
         * @param event_data Event data structure (not used in this selector)
         * @return Index of selected vertex, or -1 if none found
         */
        virtual int selectVertex(larcv::IOManager* larcv_io, 
                                larlite::storage_manager* larlite_io, 
                                gen2ntuple::EventData* event_data,
                                gen2ntuple::RecoData* reco_data) override;
    };

}

#endif