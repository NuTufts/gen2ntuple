#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

namespace gen2ntuple {

/**
 * @brief Simple logging utility
 */
class Logger {
public:
    enum Level {
        ERROR = 0,
        WARNING = 1,
        INFO = 2,
        DEBUG = 3
    };
    
    /**
     * @brief Get singleton instance
     */
    static Logger& getInstance();
    
    /**
     * @brief Set logging level
     */
    void setLevel(Level level);
    
    /**
     * @brief Set log file (optional, defaults to stdout)
     */
    void setLogFile(const std::string& filename);
    
    /**
     * @brief Log a message
     */
    void log(Level level, const std::string& message);
    
    /**
     * @brief Convenience methods
     */
    void error(const std::string& message);
    void warning(const std::string& message);
    void info(const std::string& message);
    void debug(const std::string& message);

private:
    Logger() = default;
    ~Logger() = default;
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    
    Level level_ = INFO;
    std::unique_ptr<std::ofstream> log_file_;
    
    std::ostream& getStream();
    const char* levelToString(Level level);
};

// Convenience macros
#define LOG_ERROR(msg) Logger::getInstance().error(msg)
#define LOG_WARNING(msg) Logger::getInstance().warning(msg)  
#define LOG_INFO(msg) Logger::getInstance().info(msg)
#define LOG_DEBUG(msg) Logger::getInstance().debug(msg)

} // namespace gen2ntuple