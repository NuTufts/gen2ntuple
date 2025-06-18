#include "Logger.h"
#include <iomanip>
#include <chrono>

namespace gen2ntuple {

Logger& Logger::getInstance() {
    static Logger instance;
    return instance;
}

void Logger::setLevel(Level level) {
    level_ = level;
}

void Logger::setLogFile(const std::string& filename) {
    log_file_ = std::make_unique<std::ofstream>(filename);
    if (!log_file_->is_open()) {
        std::cerr << "Warning: Could not open log file: " << filename << std::endl;
        log_file_.reset();
    }
}

void Logger::log(Level level, const std::string& message) {
    if (level > level_) {
        return;
    }
    
    // Get current time
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    
    auto& stream = getStream();
    stream << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
    stream << "." << std::setfill('0') << std::setw(3) << ms.count();
    stream << " [" << levelToString(level) << "] " << message << std::endl;
}

void Logger::error(const std::string& message) {
    log(ERROR, message);
}

void Logger::warning(const std::string& message) {
    log(WARNING, message);
}

void Logger::info(const std::string& message) {
    log(INFO, message);
}

void Logger::debug(const std::string& message) {
    log(DEBUG, message);
}

std::ostream& Logger::getStream() {
    return log_file_ && log_file_->is_open() ? *log_file_ : std::cout;
}

const char* Logger::levelToString(Level level) {
    switch (level) {
        case ERROR: return "ERROR";
        case WARNING: return "WARN";
        case INFO: return "INFO";
        case DEBUG: return "DEBUG";
        default: return "UNKNOWN";
    }
}

} // namespace gen2ntuple