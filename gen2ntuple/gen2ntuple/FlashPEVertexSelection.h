#ifndef __GEN2NTUPLE_FLASH_PE_VERTEX_SELECTION_H__
#define __GEN2NTUPLE_FLASH_PE_VERTEX_SELECTION_H__

#include "VertexSelectionBase.h"

namespace gen2ntuple {

    /**
     * @brief Select vertex based on flash prediction matching
     * 
     * This selection method picks the vertex with the lowest value of:
     * sinkhorn_divergence * |1 + fractional_pe_error|
     * 
     * This combines the shape matching (sinkhorn divergence) with the
     * total PE matching (fractional error) to find the best overall match
     * to the observed optical flash.
     */
    class FlashPEVertexSelection : public VertexSelectionBase {
    public:
        FlashPEVertexSelection() {}
        virtual ~FlashPEVertexSelection() {}
        
        /**
         * @brief Set flash predictions for vertices
         */
        virtual void setFlashPredictions(const std::vector<std::vector<float>>& flash_preds,
                                        const std::vector<float>& sinkhorn_divs) override {
            flash_predictions_ = flash_preds;
            sinkhorn_divs_ = sinkhorn_divs;
        }
        
        /**
         * @brief Select vertex with best flash prediction match
         * 
         * @param larcv_io LArCV IO manager
         * @param larlite_io LArLite storage manager
         * @param nuvtx_v Vector of neutrino vertex candidates
         * @param event_data Event data structure (contains flash prediction results)
         * @return Index of selected vertex, or -1 if none found
         */
        virtual int selectVertex(larcv::IOManager* larcv_io, 
                                larlite::storage_manager* larlite_io, 
                                gen2ntuple::EventData* event_data,
                                gen2ntuple::RecoData* reco_data) override;
                                
    private:
        /**
         * @brief Calculate the flash matching score
         * 
         * @param sinkhorn_div Sinkhorn divergence value
         * @param pred_pe Predicted total PE
         * @param obs_pe Observed total PE
         * @return Combined score (lower is better)
         */
        float calculateFlashScore(float sinkhorn_div, float pred_pe, float obs_pe) const;
        
    private:
        std::vector<std::vector<float>> flash_predictions_;
        std::vector<float> sinkhorn_divs_;
    };

}

#endif