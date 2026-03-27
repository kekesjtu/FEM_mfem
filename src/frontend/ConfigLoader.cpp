#include "fem/frontend/ConfigLoader.hpp"

#include "nlohmann/json.hpp"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem::frontend
{
using json = nlohmann::json;

template <typename T>
T ReadOrDefault(const json &node, const std::string &key, T fallback)
{
    if (!node.contains(key))
    {
        return fallback;
    }
    return node.at(key).get<T>();
}
}// namespace fem::frontend
