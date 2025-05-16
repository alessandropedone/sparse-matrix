/**
 * @file json_utility.hpp
 * @brief Utility functions for reading and writing JSON files using nlohmann::json.
 *
 * This header provides functions to read JSON data from a file and save JSON data to a file.
 * It uses the nlohmann::json library for parsing and serialization.
 *
 * @author 
 * @date 
 */

#ifndef JSON_UTILITY_HPP
#define JSON_UTILITY_HPP
#include <nlohmann/json.hpp>
#include <fstream>



namespace json_utility
{
    using json = nlohmann::json;

    /**
     * @brief Reads and parses a JSON file.
     *
     * Opens the specified file, reads its contents, and parses it into a JSON object.
     * Throws a std::runtime_error if the file cannot be opened.
     *
     * @param filename The path to the JSON file to be read.
     * @return json The parsed JSON object.
     * @throws std::runtime_error If the file cannot be opened.
     */
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


    /**
     * @brief Saves JSON data to a file with pretty formatting.
     *
     * This function writes the provided JSON object to the specified file.
     * The output is formatted with an indentation of 4 spaces for readability.
     *
     * @param filename The path to the file where the JSON data will be saved.
     * @param data The JSON object to be saved.
     *
     * @throws std::runtime_error If the file cannot be opened for writing.
     */
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