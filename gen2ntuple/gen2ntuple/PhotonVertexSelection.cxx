#include "PhotonVertexSelection.h"

#include "TF1.h"

#include "larlite/DataFormat/storage_manager.h"
#include "larlite/DataFormat/opflash.h"
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larflow/Reco/SinkhornFlashDivergence.h"


namespace gen2ntuple {

    PhotonVertexSelection::PhotonVertexSelection()
    {

        flash_predictor.set_verbosity( larcv::msg::kINFO );
        flash_predictor.setChargeToPhotonParams(100.0, 23.6e-3, 24000.0, 0.7);
        flash_predictor.setTrackConversionParams(3, 3, 0.3, 0.5);
        flash_predictor.setShowerConversionParams(3, 3);
        adc_threshold = 10.0;

        // we define two likelihood calculators, 
        // one for first for dwall>15.0 and dwall<15.0

        fL_farwall_fracerr  = new TF1( "farwall_fracerr", "gaus(0) + [3]*exp(-[4]*x)", -2.0, 100.0 );
        farwall_frac_const = 190.0;
        farwall_frac_mean  = 1.0;
        farwall_frac_sigma = 0.7;
        farwall_frac_expconst = 150.0;
        farwall_frac_lambda    = 1.0;
        farwall_frac_gaus_norm = 1.0/( farwall_frac_sigma*sqrt(2.0*3.14159)*farwall_frac_const);
        farwall_frac_exp_norm  = farwall_frac_lambda/farwall_frac_expconst;
        fL_farwall_fracerr->SetParameter(0, farwall_frac_const);
        fL_farwall_fracerr->SetParameter(1, farwall_frac_mean);
        fL_farwall_fracerr->SetParameter(2, farwall_frac_sigma);
        fL_farwall_fracerr->SetParameter(3, farwall_frac_expconst);
        fL_farwall_fracerr->SetParameter(4, farwall_frac_lambda);

        fL_farwall_sinkdiv  = new TF1( "farwall_sinkdiv", "gaus(0) + [3]*exp(-[4]*x)", 0.0, 1000.0 );
        farwall_sink_const = 96.0;
        farwall_sink_mean  = 1.4;
        farwall_sink_sigma = 0.7;
        farwall_sink_expconst = 245.0;
        farwall_sink_lambda    = 0.17;
        farwall_sink_gaus_norm = 1.0/( farwall_sink_sigma*sqrt(2.0*3.14159)*farwall_sink_const);
        farwall_sink_exp_norm  = farwall_sink_lambda/farwall_sink_expconst;
        fL_farwall_sinkdiv->SetParameter(0, farwall_sink_const);
        fL_farwall_sinkdiv->SetParameter(1, farwall_sink_mean);
        fL_farwall_sinkdiv->SetParameter(2, farwall_sink_sigma);
        fL_farwall_sinkdiv->SetParameter(3, farwall_sink_expconst);
        fL_farwall_sinkdiv->SetParameter(4, farwall_sink_lambda);

        fL_nearwall_fracerr = new TF1( "nearwall_fracerr", "gaus(0) + [3]*exp(-[4]*x)", -2.0, 100.0 );
        nearwall_frac_const = 76.0;
        nearwall_frac_mean  = 0.84;
        nearwall_frac_sigma = 0.50;
        nearwall_frac_expconst  = 40.0;
        nearwall_frac_lambda    = 1.0;
        nearwall_frac_gaus_norm = 1.0/( nearwall_frac_sigma*sqrt(2.0*3.14159)*nearwall_frac_const);
        nearwall_frac_exp_norm  = nearwall_frac_lambda/nearwall_frac_expconst;
        fL_nearwall_fracerr->SetParameter(0, nearwall_frac_const);
        fL_nearwall_fracerr->SetParameter(1, nearwall_frac_mean);
        fL_nearwall_fracerr->SetParameter(2, nearwall_frac_sigma);
        fL_nearwall_fracerr->SetParameter(3, nearwall_frac_expconst);
        fL_nearwall_fracerr->SetParameter(4, nearwall_frac_lambda);

        fL_nearwall_fracerr = new TF1( "nearwall_sinkdiv", "[0]*exp(-[1]*x)", 0.0, 1000.0 );
        nearwall_sink_expconst  = 91.0;
        nearwall_sink_lambda    = 0.08;  
        nearwall_sink_exp_norm  = nearwall_frac_lambda/nearwall_frac_expconst;
        fL_nearwall_fracerr->SetParameter(0, nearwall_sink_expconst);
        fL_nearwall_fracerr->SetParameter(1, nearwall_sink_lambda);
    }

    PhotonVertexSelection::~PhotonVertexSelection()
    {
        delete fL_farwall_fracerr;
        delete fL_farwall_sinkdiv;
        delete fL_nearwall_fracerr;
        delete fL_nearwall_sinkdiv;
    }

    int PhotonVertexSelection::selectVertex(larcv::IOManager* larcv_io, 
                                            larlite::storage_manager* larlite_io, 
                                            gen2ntuple::EventData* event_data,
                                            gen2ntuple::RecoData* reco_data) 
    {
        // we select only vertices with showers
        // we get the vertex, we calculate the flash prediction and sinkhorn divergence
        // we then use the flash prediction and sinkhorn divergence to choose the best vertex

        if (!reco_data || !reco_data->nuvtx_v || reco_data->nuvtx_v->empty()) {
            return -1;
        }
        
        const auto& nuvtx_v = *(reco_data->nuvtx_v);
        const auto& nusel_v = *(reco_data->nusel_v);

        int nvertices = nuvtx_v.size();

        // loop through the vertices and gather the info we need
        std::vector< float > totpe_pred_v(nvertices,0);
        std::vector< float > sink_v(nvertices,-1.0);
        std::vector< float > fracerr_v( nvertices, -1.0 );
        std::vector< float > dwall_v( nvertices, -999.0 );
        std::vector< float > median_pixsum_v( nvertices, 0.0 );
        std::vector< float > likelihood( nvertices, 0.0 );
        std::vector< std::vector<float> > totalpixelsum_v( nvertices );
        std::vector< std::vector<float> > leadingshower_vv( nvertices );


        // calculate the flash predictions
        float obs_total_pe = 0.0;
        std::vector<float> obs_pe_per_pmt(32, 0.0);
        try {
            auto ev_opflash = (larlite::event_opflash*)(larlite_io->get_data(larlite::data::kOpFlash, "simpleFlashBeam"));
            if (ev_opflash && ev_opflash->size() > 0) {
                auto& flash = ev_opflash->at(0);
                int nOpDets = std::min(static_cast<size_t>(32), flash.nOpDets());
                for (int ipmt = 0; ipmt < nOpDets; ipmt++) {
                    obs_pe_per_pmt[ipmt] = flash.PE(ipmt);
                    obs_total_pe += obs_pe_per_pmt[ipmt];
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not get opflash for event " << event_data->event << ": " << e.what() << std::endl;
        }
        
        // Get ADC images for flash prediction
        auto ev_adc = (larcv::EventImage2D*)larcv_io->get_data(larcv::kProductImage2D, "wire");
        const std::vector<larcv::Image2D>& adc_v = ev_adc->as_vector();
        
        for (size_t ivtx = 0; ivtx < nuvtx_v.size(); ivtx++) {
            const auto& vtx = nuvtx_v[ivtx];

            // count primaries
            int nprim_showers = 0;
            int nprim_tracks  = 0;
            for (int ishower=0; ishower<(int)vtx.shower_v.size(); ishower++) {
                if ( vtx.shower_isSecondary_v.at(ishower)==0 )
                    nprim_showers++;
            }
            for (int itrack=0; itrack<(int)vtx.track_v.size(); itrack++) {
                if ( vtx.track_isSecondary_v.at(itrack)==0 ){
                    nprim_tracks++;
                }
            }

            if ( nprim_showers==0 )
                continue;

            // determine dwall variable values and closest vertices for these metrics
            std::vector<float> reco_vtx = {vtx.pos[0], vtx.pos[1], vtx.pos[2]};
            float x_dwall_reco_nuvtx = dwall( reco_vtx[0], reco_vtx[1], reco_vtx[2] );
            dwall_v[ivtx] = x_dwall_reco_nuvtx;
            
            // Calculate flash prediction variables
            try {
                auto predicted_flash = flash_predictor.predictFlash(
                    vtx,
                    adc_v,
                    adc_threshold,
                    true,  // use trilinear interpolation
                    false  // include all particles
                );
                
                float pred_total_pe = flash_predictor.getTotalPredictedPE();
                auto pe_per_pmt_map = flash_predictor.getPredictedPE();
                
                // Convert map to vector for Sinkhorn calculation
                std::vector<float> pred_pe_per_pmt(32, 0.0);
                float x_predicted_totpe = 0.;
                for (int ipmt = 0; ipmt < 32; ipmt++) {
                    if (pe_per_pmt_map.find(ipmt) != pe_per_pmt_map.end()) {
                        pred_pe_per_pmt[ipmt] = pe_per_pmt_map[ipmt];
                        x_predicted_totpe += pred_pe_per_pmt[ipmt];
                    }
                }
                totpe_pred_v[ivtx] = x_predicted_totpe;
                
                // Calculate flash prediction metrics
                float x_sinkhorndiv  = calculateSinkhornDivergence(pred_pe_per_pmt, obs_pe_per_pmt);
                float x_totpefracerr = (pred_total_pe - obs_total_pe) / (0.1 + obs_total_pe);
                sink_v[ivtx]    = x_sinkhorndiv;
                fracerr_v[ivtx] = x_totpefracerr;

                float Lsink = 0.;
                float Lfrac = 0.;
                if ( x_dwall_reco_nuvtx>20.0 ) {
                    Lfrac = fL_farwall_fracerr->Eval( 1.0+x_totpefracerr );
                    Lsink = fL_farwall_sinkdiv->Eval( x_sinkhorndiv );
                }
                else {
                    Lfrac = fL_nearwall_fracerr->Eval( 1.0+x_totpefracerr );
                    Lsink = fL_nearwall_sinkdiv->Eval( x_sinkhorndiv );
                }

                float x_likelihood = Lfrac*Lsink;
                likelihood[ivtx] = x_likelihood;

                // Calculate total pixel sum across planes
		        std::vector<float> leadingshower_v(3,-1.0);
                std::vector<float> totalpixelsum = calculateProngPixelSum( *larcv_io, vtx, leadingshower_v, true );
                totalpixelsum_v[ivtx]  = totalpixelsum;
		        leadingshower_vv[ivtx] = leadingshower_v;
                std::sort( totalpixelsum.begin(), totalpixelsum.end() );
                float x_allprong_median_pixsum = totalpixelsum[1];
		        std::sort( leadingshower_v.begin(), leadingshower_v.end() );
		        float x_reco_median_pixsum = leadingshower_v[1];
                median_pixsum_v[ivtx] = x_reco_median_pixsum;
                
            } catch (const std::exception& e) {
                std::cerr << "Warning: Flash prediction failed for entry " << event_data->event << ", "
                          << " vertex " << ivtx << ": " << e.what() << std::endl;
            }
        
        }//end of vertex loop

        // decide

        return 0;
    }

    float PhotonVertexSelection::dwall( float x, float y, float z )
    {

        float xmin = x;
        float xmax = 256.0-x;
        float ymin = y + 116.5; // y - -116.5
        float ymax = 116.5-y;
        float zmin = z;
        float zmax = 1036.0-z;

        std::array<float,6> dists = { xmin, xmax, ymin, ymax, zmin, zmax };

        if ( xmin>=0 && xmax>=0 && ymin>=0 && ymax>=0 && zmin>=0 && zmax>=0 ) {
            // is inside, find min
            float mindist = 1.0e9;
            for (auto& dist : dists ) {
                if ( dist < mindist )
                    mindist = dist;
            }
            return mindist;
        }
        else {
            // outside: only negative values
            float mindist = 1.0e9;
            for ( auto& dist : dists ) {
                if ( dist>=0 ) continue; // must be on the "wrong" side of the boundary
                if ( fabs(dist)<mindist )
                    mindist = fabs(dist);
            }
            // return a negative value to indicate we are outside the TPC
            return -mindist;
        }

        // should never get here
        return 1.0e9;
    }

    float PhotonVertexSelection::calculateSinkhornDivergence( const std::vector<float>& pred_pe, 
                                                              const std::vector<float>& obs_pe) {
        try {
            larflow::reco::SinkhornFlashDivergence sinkhorn_calc;
            return sinkhorn_calc.calculateDivergence(pred_pe, obs_pe, 1.0, 100, 1e-6);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Sinkhorn calculation failed: " << e.what() << std::endl;
            return -999.0;
        }
        return -999.0;
    }

    std::vector<float> PhotonVertexSelection::calculateProngPixelSum( 
        larcv::IOManager& larcv_io,
        const larflow::reco::NuVertexCandidate& nuvtx,
        std::vector<float>& leading_shower_pixsum_v,
        bool primary_only ) 
    {

        int ntracks  = nuvtx.track_v.size();
        int nshowers = nuvtx.shower_v.size();

        std::vector<float> total_charge(3,0.0);
        leading_shower_pixsum_v.resize(3,-1.0);
        for (int p=0; p<3; p++)
          leading_shower_pixsum_v[p] = -1.0;

        auto ev_img = (larcv::EventImage2D*)larcv_io.get_data("image2d", "wire");
        auto image2Dvec = ev_img->as_vector();
        float threshold = 10.0f; // ADC threshold

        for (int p = 0; p < 3; ++p ) {

            auto const& img = image2Dvec.at(p);
            std::set< std::array<int,3> > pixels_visited; // prevent double counting pixels

            for (int track_idx=0; track_idx<ntracks; track_idx++) {

                float clusterCharge = 0.0f;

                if ( primary_only && nuvtx.track_isSecondary_v[track_idx]==1 ) {
                    // is secondary and we are only including primaries
                    continue;
                }

                auto const& cluster = nuvtx.track_hitcluster_v.at(track_idx);

                for (auto const& hit : cluster) {

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

                total_charge[p] += clusterCharge;

            }//end of track loop

            // loop over showers and get pixelsums for the plane
       	    int leading_shower_index = -1;
       	    float leading_shower_pixsum = -1.0;
            for (int shower_idx=0; shower_idx<nshowers; shower_idx++ ) {

                float clusterCharge = 0.0f;

                if ( primary_only && nuvtx.shower_isSecondary_v[shower_idx]==1 ) {
                    continue;
                }

                auto const& cluster = nuvtx.shower_v.at(shower_idx);

                for (auto const& hit : cluster) {
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
                }// loop over hits

    	    if ( leading_shower_pixsum<0 || clusterCharge>leading_shower_pixsum ) {
    	      leading_shower_pixsum = clusterCharge;
    	      leading_shower_index = shower_idx;
    	    }

                total_charge[p] += clusterCharge;

            }// loop over showers

    	leading_shower_pixsum_v[p] = leading_shower_pixsum;
    
        }// loop over planes

        return total_charge;
    }

}