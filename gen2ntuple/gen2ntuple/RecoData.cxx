#include "RecoData.h"

namespace gen2ntuple {

    void RecoData::setBranchAddresses( TChain* recotree )
    {

        recotree->SetBranchAddress("run",        &run);
        recotree->SetBranchAddress("subrun",     &subrun);
        recotree->SetBranchAddress("event",      &event);
        
        recotree->SetBranchAddress("nuvetoed_v", &nuvtx_v );
        recotree->SetBranchAddress("nu_sel_v",   &nusel_v);
    }


}