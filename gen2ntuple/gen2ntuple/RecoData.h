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

        RecoData() {};
        ~RecoData() {};

        int run;
        int subrun;
        int event;
        std::vector<larflow::reco::NuVertexCandidate>*    nuvtx_v;
        std::vector<larflow::reco::NuSelectionVariables>* nusel_v;

        void setBranchAddresses( TChain* recotree );


    };

}

#endif