#pragma once

#include "EventData.h"
#include "Configuration.h"

namespace gen2ntuple {

/**
 * @brief Main event processing coordinator
 * 
 * This class orchestrates the processing of individual events,
 * coordinating between the various processing modules.
 */
class EventProcessor {
public:
    /**
     * @brief Constructor
     */
    EventProcessor(const Configuration& config);
    
    /**
     * @brief Process a single event
     * @param event_data Event data structure to fill
     * @return true if processing successful
     */
    bool processEvent(EventData& event_data);

private:
    const Configuration& config_;
    
    // TODO: Add processing module instances
};

} // namespace gen2ntuple