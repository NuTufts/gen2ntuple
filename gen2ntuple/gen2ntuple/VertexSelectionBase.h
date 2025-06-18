#ifndef __GEN2NTUPLE_VERTEX_SELECTION_BASE_H__
#define __GEN2NTUPLE_VERTEX_SELECTION_BASE_H__

#include "larlite/DataFormat/storage_manager.h"
#include "larcv/core/DataFormat/IOManager.h"
#include "larflow/Reco/NuVertexCandidate.h"
#include "EventData.h"
#include "RecoData.h"

namespace gen2ntuple {

    class VertexSelectionBase {
        public:

        VertexSelectionBase() {};
        virtual ~VertexSelectionBase() {};
        
        /**
         * @brief Set flash predictions for vertices (used by flash-based selectors)
         */
        virtual void setFlashPredictions(const std::vector<std::vector<float>>& flash_preds,
                                        const std::vector<float>& sinkhorn_divs) {
            // Default implementation does nothing
        }

        /**
         * @brief Method that child classes must implement
         * 
         * Should return the index of the NuVertexCandidate container that
         * is accepted as the neutrino interaction for the ntuple.
         */
        virtual int selectVertex( larcv::IOManager* larcv_io, 
                          larlite::storage_manager* larlite_io, 
                          gen2ntuple::EventData* event_data,
                          gen2ntuple::RecoData* reco_data ) = 0;

    };

}

#endif