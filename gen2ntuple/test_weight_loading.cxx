// Test program for MCTruthProcessor weight loading functionality

#include <iostream>
#include <string>
#include "gen2ntuple/MCTruthProcessor.h"
#include "gen2ntuple/EventData.h"

int main(int argc, char* argv[]) {
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <weight_file.pkl>" << std::endl;
        return 1;
    }
    
    std::string weight_file = argv[1];
    
    // Create MCTruthProcessor instance
    gen2ntuple::MCTruthProcessor mc_processor;
    
    // Set weight file
    mc_processor.setWeightFile(weight_file);
    
    // Load weights
    std::cout << "Loading weights from: " << weight_file << std::endl;
    if (!mc_processor.loadWeights()) {
        std::cerr << "Failed to load weights!" << std::endl;
        return 1;
    }
    
    // Test a few weight lookups
    // These are example run/subrun/event numbers - you might need to adjust
    // based on what's actually in your weight files
    int test_cases[][3] = {
        {5226, 1, 50},    // Example from run 1
        {5226, 10, 500},  // Another example
        {9999, 1, 1}      // Non-existent event (should return 1.0)
    };
    
    for (int i = 0; i < 3; i++) {
        int run = test_cases[i][0];
        int subrun = test_cases[i][1];
        int event = test_cases[i][2];
        
        float weight = mc_processor.getEventWeight(run, subrun, event);
        std::cout << "Weight for run=" << run 
                  << " subrun=" << subrun 
                  << " event=" << event 
                  << " : " << weight << std::endl;
    }
    
    std::cout << "Weight loading test completed successfully!" << std::endl;
    
    return 0;
}