#pragma once

#include "EventData.h"
#include "TTree.h"
#include "TFile.h"
#include <memory>

namespace gen2ntuple {

/**
 * @brief Manages ROOT tree creation and branch setup
 * 
 * This class handles the creation of ROOT trees and sets up all branches
 * to point to the appropriate members of EventData and POTData structures.
 * It provides a clean interface for tree management and ensures consistent
 * branch naming and types.
 */
class BranchManager {
public:
    /**
     * @brief Constructor
     * @param output_file Output ROOT file
     * @param is_mc Whether this is MC data (affects which branches are created)
     * @param include_keypoints Whether to include keypoint branches
     */
    BranchManager(TFile* output_file, bool is_mc = false, bool include_keypoints = false);
    
    /**
     * @brief Destructor
     */
    ~BranchManager();
    
    /**
     * @brief Set up all branches for the event tree
     * @param event_data Pointer to EventData structure
     */
    void setupEventBranches(EventData* event_data);
    
    /**
     * @brief Set up branches for the POT tree (MC only)
     * @param pot_data Pointer to POTData structure
     */
    void setupPOTBranches(POTData* pot_data);
    
    /**
     * @brief Fill the event tree
     */
    void fillEventTree();
    
    /**
     * @brief Fill the POT tree
     */
    void fillPOTTree();
    
    /**
     * @brief Write trees to file
     */
    void write();
    
    /**
     * @brief Get the event tree
     */
    TTree* getEventTree() { return event_tree_; }
    
    /**
     * @brief Get the POT tree
     */
    TTree* getPOTTree() { return pot_tree_; }
    
    /**
     * @brief Get number of entries in event tree
     */
    Long64_t getEventEntries() const;
    
    /**
     * @brief Get number of entries in POT tree
     */
    Long64_t getPOTEntries() const;

private:
    TFile* output_file_;
    bool is_mc_;
    bool include_keypoints_;
    
    TTree* event_tree_;
    TTree* pot_tree_;
    
    /**
     * @brief Set up event identification branches
     */
    void setupEventIDBranches(EventData* data);
    
    /**
     * @brief Set up MC truth branches
     */
    void setupMCTruthBranches(EventData* data);
    
    /**
     * @brief Set up vertex reconstruction branches
     */
    void setupVertexBranches(EventData* data);
    
    /**
     * @brief Set up keypoint branches
     */
    void setupKeypointBranches(EventData* data);
    
    /**
     * @brief Set up track branches
     */
    void setupTrackBranches(EventData* data);
    
    /**
     * @brief Set up shower branches
     */
    void setupShowerBranches(EventData* data);
    
    /**
     * @brief Set up event-level feature branches
     */
    void setupEventFeatureBranches(EventData* data);
    
    /**
     * @brief Helper function to create array branch
     * @param tree Tree to add branch to
     * @param name Branch name
     * @param address Data address
     * @param size_var Variable containing array size
     * @param max_size Maximum array size
     * @param type Data type character (I, F, O, etc.)
     */
    template<typename T>
    void createArrayBranch(TTree* tree, const char* name, T* address, 
                          const char* size_var, int max_size, const char* type);
};

} // namespace gen2ntuple