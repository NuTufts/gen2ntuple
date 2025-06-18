#pragma once

#include <string>
#include <vector>

namespace gen2ntuple {

/**
 * @brief Configuration class for runtime parameters
 * 
 * This class handles command line parsing and provides access to all
 * runtime configuration parameters for the ntuple creation process.
 */
class Configuration {
public:
    /**
     * @brief Constructor
     */
    Configuration();
    
    /**
     * @brief Parse command line arguments
     * @param argc Number of arguments
     * @param argv Argument array
     * @return true if parsing successful, false otherwise
     */
    bool parseCommandLine(int argc, char** argv);
    
    /**
     * @brief Print usage information
     */
    void printUsage(const char* program_name) const;
    
    /**
     * @brief Validate configuration
     * @return true if configuration is valid
     */
    bool validate() const;
    
    // ===============================
    // Input/Output Configuration
    // ===============================
    
    /**
     * @brief Get input KPSReco files
     */
    const std::vector<std::string>& getKPSRecoFiles() const { return kpsreco_files_; }
    
    /**
     * @brief Get truth file list or single truth file
     */
    const std::string& getTruthFile() const { return truth_file_; }
    
    /**
     * @brief Get output filename
     */
    const std::string& getOutputFile() const { return output_file_; }
    
    /**
     * @brief Get weight file (for MC)
     */
    const std::string& getWeightFile() const { return weight_file_; }
    
    /**
     * @brief Get CNN model path
     */
    const std::string& getModelPath() const { return model_path_; }
    
    // ===============================
    // Processing Configuration
    // ===============================
    
    /**
     * @brief Check if this is MC data
     */
    bool isMC() const { return is_mc_; }
    
    /**
     * @brief Check if using dlana input files
     */
    bool isDLAna() const { return is_dlana_; }
    
    /**
     * @brief Check if keypoint processing is disabled
     */
    bool isKeypointsDisabled() const { return no_keypoints_; }
    
    /**
     * @brief Get device for CNN processing
     */
    const std::string& getDevice() const { return device_; }
    
    /**
     * @brief Get maximum number of events to process
     */
    int getMaxEvents() const { return max_events_; }
    
    /**
     * @brief Get starting event number
     */
    int getStartEvent() const { return start_event_; }
    
    /**
     * @brief Get number of worker threads
     */
    int getNumThreads() const { return num_threads_; }
    
    /**
     * @brief Check if verbose output is enabled
     */
    bool isVerbose() const { return verbose_; }
    
    /**
     * @brief Check if debug output is enabled
     */
    bool isDebug() const { return debug_; }

private:
    // Input/Output
    std::vector<std::string> kpsreco_files_;
    std::string truth_file_;
    std::string output_file_;
    std::string weight_file_;
    std::string model_path_;
    
    // Processing flags
    bool is_mc_;
    bool is_dlana_;
    bool no_keypoints_;
    std::string device_;
    
    // Event processing
    int max_events_;
    int start_event_;
    int num_threads_;
    
    // Output control
    bool verbose_;
    bool debug_;
    
    /**
     * @brief Set default values
     */
    void setDefaults();
};

} // namespace gen2ntuple