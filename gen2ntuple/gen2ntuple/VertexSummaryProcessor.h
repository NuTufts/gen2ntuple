#ifndef __GEN2NTUPLE_VERTEX_SUMMARY_PROCESSOR_H__
#define __GEN2NTUPLE_VERTEX_SUMMARY_PROCESSOR_H__

#include "EventData.h"
#include "RecoData.h"

namespace gen2ntuple {

    class VertexSummaryProcessor {

    public:

        VertexSummaryProcessor() {};
        ~VertexSummaryProcessor() {};


        bool calculateHitAndChargeFractions( EventData* event_data );

        bool calculateEventPCA( EventData* event_data, RecoData* reco_data );

    };


}

#endif