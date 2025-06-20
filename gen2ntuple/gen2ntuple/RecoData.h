#ifndef __GEN2NTUPLE_RECODATA_H__
#define __GEN2NTUPLE_RECODATA_H__


#include <vector>
#include "TChain.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/NuSelectionVariables.h"
#include "larflow/Reco/KPCluster.h"

/**
 * @brief Interface to data from the KPSRecoManagerTree
 * 
 */



namespace gen2ntuple {

    class RecoData {

    public:

        RecoData() 
        : nuvtx_v(nullptr), 
        nusel_v(nullptr),
        kpc_nu_v(nullptr),
        kpc_track_v(nullptr),
        kpc_shower_v(nullptr),
        kpc_cosmic_v(nullptr) 
        {};
        ~RecoData() {
            if (nuvtx_v)  delete nuvtx_v;
            if (nusel_v)  delete nusel_v;
            if (kpc_nu_v) delete kpc_nu_v;
            if (kpc_track_v)  delete kpc_track_v;
            if (kpc_shower_v) delete kpc_shower_v;
            if (kpc_cosmic_v) delete kpc_cosmic_v;
        };

        // int run;
        // int subrun;
        // int event;
        std::vector<larflow::reco::NuVertexCandidate>*    nuvtx_v;
        std::vector<larflow::reco::NuSelectionVariables>* nusel_v;
        std::vector<larflow::reco::KPCluster>*            kpc_nu_v;
        std::vector<larflow::reco::KPCluster>*            kpc_track_v;
        std::vector<larflow::reco::KPCluster>*            kpc_shower_v;
        std::vector<larflow::reco::KPCluster>*            kpc_cosmic_v;

        void setBranchAddresses( TChain* recotree );


    };

}

#endif