#ifndef __GEN2NTUPLE_PHOTON_VERTEX_SELECTION_H__
#define __GEN2NTUPLE_PHOTON_VERTEX_SELECTION_H__

#include "VertexSelectionBase.h"

#include "larflow/Reco/NuVertexFlashPrediction.h"

class TF1;

namespace gen2ntuple {

    /**
     * @brief Select vertex based on keypoint type priority and visible energy
     * 
     * This selection method prioritizes neutrino keypoints (type==0).
     * If multiple nu keypoints found, picks the one with highest visible energy.
     * If no nu keypoints found, picks non-nu vertex with highest visible energy.
     * Corresponds to the Python highest_kprank_with_visenergy() function.
     */
    class PhotonVertexSelection : public VertexSelectionBase {
    public:
        PhotonVertexSelection();
        virtual ~PhotonVertexSelection();
        
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
                                
    protected:

        float adc_threshold;
        TF1* fL_farwall_fracerr;
        TF1* fL_farwall_sinkdiv;
        TF1* fL_nearwall_fracerr;
        TF1* fL_nearwall_sinkdiv;

        float farwall_frac_const;
        float farwall_frac_mean;
        float farwall_frac_sigma;
        float farwall_frac_expconst;
        float farwall_frac_lambda;
        float farwall_frac_gaus_norm;
        float farwall_frac_exp_norm;

        float farwall_sink_const;
        float farwall_sink_mean;
        float farwall_sink_sigma;
        float farwall_sink_expconst;
        float farwall_sink_lambda;
        float farwall_sink_gaus_norm;
        float farwall_sink_exp_norm;

        float nearwall_frac_const;
        float nearwall_frac_mean;
        float nearwall_frac_sigma;
        float nearwall_frac_expconst;
        float nearwall_frac_lambda;
        float nearwall_frac_gaus_norm;
        float nearwall_frac_exp_norm;

        float nearwall_sink_expconst;
        float nearwall_sink_lambda; 
        float nearwall_sink_exp_norm;

        larflow::reco::NuVertexFlashPrediction flash_predictor;

        float dwall( float x, float y, float z );
        float calculateSinkhornDivergence( const std::vector<float>& pred_pe, 
                                           const std::vector<float>& obs_pe);
        std::vector<float> calculateProngPixelSum( 
            larcv::IOManager& larcv_io,
            const larflow::reco::NuVertexCandidate& nuvtx,
            std::vector<float>& leading_shower_pixsum_v,
            bool primary_only );

    };

}

#endif