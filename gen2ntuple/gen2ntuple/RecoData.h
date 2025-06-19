#ifndef __GEN2NTUPLE_RECODATA_H__
#define __GEN2NTUPLE_RECODATA_H__


#include <vector>
#include "TChain.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "larflow/Reco/NuSelectionVariables.h"

/**
 * @brief Interface to data from the KPSRecoManagerTree
 * 
 */



namespace gen2ntuple {

    class RecoData {

    public:

        RecoData() : nuvtx_v(nullptr), nusel_v(nullptr) {};
        ~RecoData() {
            if (nuvtx_v) delete nuvtx_v;
            if (nusel_v) delete nusel_v;
        };

        // int run;
        // int subrun;
        // int event;
        std::vector<larflow::reco::NuVertexCandidate>*    nuvtx_v;
        std::vector<larflow::reco::NuSelectionVariables>* nusel_v;

        void setBranchAddresses( TChain* recotree );


    };

}

#endif