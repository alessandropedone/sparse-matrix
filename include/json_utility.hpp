#ifndef JSON_UTILITY_HPP
#define JSON_UTILITY_HPP
#include <nlohmann/json.hpp>
#include <fstream>

namespace json_utility
{
    using json = nlohmann::json;

    json read_json(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filename);
        }
        json data = json::parse(file);
        return data;
    }
    void save_json(const std::string &filename, const json &data)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filename);
        }
        file << data.dump(4);
        file.close();
    }
}

#endif // JSON_UTILITY_HPP