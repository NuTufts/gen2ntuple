#include "Configuration.h"
#include <iostream>
#include <getopt.h>
#include <sstream>
#include <sys/stat.h>

namespace gen2ntuple {

Configuration::Configuration() {
    setDefaults();
}

void Configuration::setDefaults() {
    kpsreco_files_.clear();
    truth_file_ = "";
    output_file_ = "dlgen2_flat_ntuple.root";
    weight_file_ = "none";
    model_path_ = "";
    vertex_selection_ = "";
    
    is_mc_ = false;
    is_dlana_ = false;
    no_keypoints_ = false;
    device_ = "cpu";
    
    max_events_ = -1;  // Process all events
    start_event_ = 0;
    num_threads_ = 1;
    
    verbose_ = false;
    debug_ = false;
}

bool Configuration::parseCommandLine(int argc, char** argv) {
    
    // Define long options
    static struct option long_options[] = {
        {"files",      required_argument, 0, 'f'},
        {"truth",      required_argument, 0, 't'},
        {"output",     required_argument, 0, 'o'},
        {"weightfile", required_argument, 0, 'w'},
        {"model",      required_argument, 0, 'm'},
        {"device",     required_argument, 0, 'd'},
        {"max-events", required_argument, 0, 'n'},
        {"start",      required_argument, 0, 's'},
        {"threads",    required_argument, 0, 'j'},
        {"vertex-selection", required_argument, 0, 'x'},
        {"mc",         no_argument,       0, 1001},
        {"dlana",      no_argument,       0, 1002},
        {"no-kp",      no_argument,       0, 1003},
        {"verbose",    no_argument,       0, 'v'},
        {"debug",      no_argument,       0, 1004},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "f:t:o:w:m:d:n:s:j:x:vh", long_options, &option_index)) != -1) {
        switch (c) {
            case 'f':
                // Parse space-separated list of files
                {
                    std::string files_str(optarg);
                    std::istringstream iss(files_str);
                    std::string file;
                    while (iss >> file) {
                        kpsreco_files_.push_back(file);
                    }
                }
                break;
                
            case 't':
                truth_file_ = optarg;
                break;
                
            case 'o':
                output_file_ = optarg;
                break;
                
            case 'w':
                weight_file_ = optarg;
                break;
                
            case 'm':
                model_path_ = optarg;
                break;
                
            case 'd':
                device_ = optarg;
                break;
                
            case 'n':
                max_events_ = std::atoi(optarg);
                break;
                
            case 's':
                start_event_ = std::atoi(optarg);
                break;
                
            case 'j':
                num_threads_ = std::atoi(optarg);
                break;

            case 'x':
                vertex_selection_ = std::string(optarg);
                
            case 1001:  // --mc
                is_mc_ = true;
                break;
                
            case 1002:  // --dlana
                is_dlana_ = true;
                break;
                
            case 1003:  // --no-kp
                no_keypoints_ = true;
                break;
                
            case 'v':
                verbose_ = true;
                break;
                
            case 1004:  // --debug
                debug_ = true;
                verbose_ = true;  // Debug implies verbose
                break;
                
            case 'h':
                printUsage(argv[0]);
                return false;
                
            case '?':
                std::cerr << "Unknown option or missing argument" << std::endl;
                printUsage(argv[0]);
                return false;
                
            default:
                std::cerr << "Unexpected option: " << c << std::endl;
                return false;
        }
    }
    
    // Handle additional positional arguments for files if -f wasn't used
    if (kpsreco_files_.empty() && optind < argc) {
        for (int i = optind; i < argc; i++) {
            kpsreco_files_.push_back(argv[i]);
        }
    }
    
    return validate();
}

void Configuration::printUsage(const char* program_name) const {
    std::cout << "Make Flat NTuples for DLGen2 Analyses\n\n";
    std::cout << "Usage: " << program_name << " [options]\n\n";
    
    std::cout << "Required arguments:\n";
    std::cout << "  -f, --files FILE...        Input KPSReco files (space-separated)\n";
    std::cout << "  -t, --truth FILE           Text file with merged_dlreco list or single merged_dlreco file\n";
    std::cout << "  -m, --model PATH           Path to prong CNN checkpoint file\n\n";
    
    std::cout << "Optional arguments:\n";
    std::cout << "  -o, --output FILE          Output file name (default: dlgen2_flat_ntuple.root)\n";
    std::cout << "  -w, --weightfile FILE      Weights file for MC (pickled python dict, default: none)\n";
    std::cout << "  -d, --device DEVICE        Device for CNN processing (cpu/cuda, default: cpu)\n";
    std::cout << "  -n, --max-events N         Maximum number of events to process (default: all)\n";
    std::cout << "  -s, --start N              Starting event number (default: 0)\n";
    std::cout << "  -j, --threads N            Number of worker threads (default: 1)\n\n";
    
    std::cout << "Data type flags:\n";
    std::cout << "  --mc                       Running over MC input\n";
    std::cout << "  --dlana                    Using merged_dlana input files\n\n";
    
    std::cout << "Processing flags:\n";
    std::cout << "  --no-kp                    Disable keypoint processing\n";
    std::cout << "  -v, --verbose              Enable verbose output\n";
    std::cout << "  --debug                    Enable debug output (implies verbose)\n";
    std::cout << "  -h, --help                 Show this help message\n\n";
    
    std::cout << "Examples:\n";
    std::cout << "  # Process MC data with CNN\n";
    std::cout << "  " << program_name << " --mc \\\n";
    std::cout << "    -f kpsreco_file1.root kpsreco_file2.root \\\n";
    std::cout << "    -t merged_dlreco_list.txt \\\n";
    std::cout << "    -m /path/to/model.pt \\\n";
    std::cout << "    -w weights.pkl \\\n";
    std::cout << "    -o output.root\n\n";
    
    std::cout << "  # Process data without keypoints\n";
    std::cout << "  " << program_name << " --no-kp \\\n";
    std::cout << "    -f kpsreco_file.root \\\n";
    std::cout << "    -t merged_dlreco_file.root \\\n";
    std::cout << "    -m /path/to/model.pt\n";
}

bool Configuration::validate() const {
    bool valid = true;
    
    // Check required arguments
    if (kpsreco_files_.empty()) {
        std::cerr << "Error: No input KPSReco files specified" << std::endl;
        valid = false;
    }
    
    if (truth_file_.empty()) {
        std::cerr << "Error: No truth file specified" << std::endl;
        valid = false;
    }
    
    if (model_path_.empty()) {
        std::cerr << "Error: No CNN model path specified" << std::endl;
        valid = false;
    }
    
    // Check file existence using stat
    auto file_exists = [](const std::string& path) {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    };
    
    for (const auto& file : kpsreco_files_) {
        if (!file_exists(file)) {
            std::cerr << "Error: KPSReco file does not exist: " << file << std::endl;
            valid = false;
        }
    }
    
    if (!file_exists(truth_file_)) {
        std::cerr << "Error: Truth file does not exist: " << truth_file_ << std::endl;
        valid = false;
    }
    
    if (!file_exists(model_path_)) {
        std::cerr << "Error: Model file does not exist: " << model_path_ << std::endl;
        valid = false;
    }
    
    // Check weight file if MC
    if (is_mc_ && weight_file_ != "none" && !file_exists(weight_file_)) {
        std::cerr << "Error: Weight file does not exist: " << weight_file_ << std::endl;
        valid = false;
    }
    
    // Check numeric parameters
    if (max_events_ == 0) {
        std::cerr << "Error: max-events cannot be zero" << std::endl;
        valid = false;
    }
    
    if (start_event_ < 0) {
        std::cerr << "Error: start-event cannot be negative" << std::endl;
        valid = false;
    }
    
    if (num_threads_ <= 0) {
        std::cerr << "Error: threads must be positive" << std::endl;
        valid = false;
    }
    
    // Check device
    if (device_ != "cpu" && device_ != "cuda") {
        std::cerr << "Error: device must be 'cpu' or 'cuda'" << std::endl;
        valid = false;
    }
    
    return valid;
}

} // namespace gen2ntuple