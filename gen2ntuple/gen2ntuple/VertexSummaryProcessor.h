#ifndef __GEN2NTUPLE_VERTEX_SUMMARY_PROCESSOR_H__
#define __GEN2NTUPLE_VERTEX_SUMMARY_PROCESSOR_H__

#include "EventData.h"

namespace gen2ntuple {

    class VertexSummaryProcessor {

    public:

        VertexSummaryProcessor() {};
        ~VertexSummaryProcessor() {};


        bool calculateHitAndChargeFractions( EventData* event_data );


    };


}

#endif