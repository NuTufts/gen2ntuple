#ifndef __GEN2NTUPLE_GET_NUCANDIDATE_INTIME_CHARGE_H__
#define __GEN2NTUPLE_GET_NUCANDIDATE_INTIME_CHARGE_H__

#include "VertexSelectionBase.h"

namespace gen2ntuple {

    /**
     * @brief Select vertex based on highest in-time reco fraction
     * 
     * This selection method prioritizes vertices based on keypoint type order
     * (0, 3, 1, 2) and selects the one with the lowest unreconstructed fraction
     * (highest reconstructed fraction) of in-time charge.
     * Corresponds to the Python highest_intime_reco_frac() function.
     */
    class GetNuCandidateIntimeCharge : public VertexSelectionBase {
    public:
        GetNuCandidateIntimeCharge() : prioritize_by_keypoint_(true) {}
        virtual ~GetNuCandidateIntimeCharge() {}
        
        /**
         * @brief Set whether to prioritize by keypoint type
         */
        void setPrioritizeByKeypoint(bool prioritize) { prioritize_by_keypoint_ = prioritize; }
        
        /**
         * @brief Select vertex with highest in-time reconstructed charge fraction
         * 
         * @param larcv_io LArCV IO manager  
         * @param larlite_io LArLite storage manager
         * @param nuvtx_v Vector of neutrino vertex candidates
         * @param event_data Event data structure (used to access nuselvar_v)
         * @return Index of selected vertex, or -1 if none found
         */
        virtual int selectVertex(larcv::IOManager* larcv_io, 
                                larlite::storage_manager* larlite_io, 
                                gen2ntuple::EventData* event_data,
                                gen2ntuple::RecoData* reco_data) override;
                                
    private:
        bool prioritize_by_keypoint_;
        
        /**
         * @brief Calculate in-time charge for a vertex candidate
         * 
         * This method masks the tracks/showers from the vertex and computes
         * the charge that is not tagged as cosmic.
         */
        float calculateIntimeCharge(const larflow::reco::NuVertexCandidate& vtx,
                                   larcv::IOManager* larcv_io) const;
    };

}

#endif