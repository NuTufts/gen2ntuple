#include "RecoData.h"

namespace gen2ntuple {

    void RecoData::setBranchAddresses( TChain* recotree )
    {
        // // Initialize the vector pointers
        // nuvtx_v = new std::vector<larflow::reco::NuVertexCandidate>();
        // nusel_v = new std::vector<larflow::reco::NuSelectionVariables>();

        // recotree->SetBranchAddress("run",        &run);
        // recotree->SetBranchAddress("subrun",     &subrun);
        // recotree->SetBranchAddress("event",      &event);
        
        recotree->SetBranchAddress("nuvetoed_v",    &nuvtx_v );
        recotree->SetBranchAddress("nu_sel_v",      &nusel_v);
        recotree->SetBranchAddress("kpc_nu_v",      &kpc_nu_v );
        recotree->SetBranchAddress("kpc_track_v",   &kpc_track_v );
        recotree->SetBranchAddress("kpc_shower_v",  &kpc_shower_v );
        recotree->SetBranchAddress("kpc_cosmic_v",  &kpc_cosmic_v );

    }


}