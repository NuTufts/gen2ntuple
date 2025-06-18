#ifndef __GEN2NTUPLE_HIGHEST_KPRANK_WITH_VISENERGY_H__
#define __GEN2NTUPLE_HIGHEST_KPRANK_WITH_VISENERGY_H__

#include "VertexSelectionBase.h"

namespace gen2ntuple {

    /**
     * @brief Select vertex based on keypoint type priority and visible energy
     * 
     * This selection method prioritizes neutrino keypoints (type==0).
     * If multiple nu keypoints found, picks the one with highest visible energy.
     * If no nu keypoints found, picks non-nu vertex with highest visible energy.
     * Corresponds to the Python highest_kprank_with_visenergy() function.
     */
    class HighestKPRankWithVisEnergy : public VertexSelectionBase {
    public:
        HighestKPRankWithVisEnergy() 
        : min_num_showers_(0), 
        min_num_tracks_(0), 
        primary_only_(true) 
        {};
        virtual ~HighestKPRankWithVisEnergy() {}
        
        /**
         * @brief Set minimum shower requirement
         */
        void setMinNumShowers(int min_showers) { min_num_showers_ = min_showers; }
        
        /**
         * @brief Set minimum track requirement
         */
        void setMinNumTracks(int min_tracks) { min_num_tracks_ = min_tracks; }

        /**
         * @brief Set if we use only primary particles to pick vertices
         */
        void setPrimaryOnly(bool doit) { primary_only_ = doit; } 
        
        /**
         * @brief Select vertex with keypoint priority and visible energy ranking
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
        int min_num_showers_;
        int min_num_tracks_;
        bool primary_only_;
    };

}

#endif