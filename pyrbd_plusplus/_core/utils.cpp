#include <pyrbd_plusplus/utils.hpp>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_set>

namespace pyrbdpp::utils
{
   
    bool isSubSet(const Set &set1, const Set &set2)
    {   
        std::unordered_set<NodeID> set2Lookup(set2.begin(), set2.end());

        // Check if set1 is a subset of set2
        for (const auto &elem : set1)
        {
            if (set2Lookup.find(elem) == set2Lookup.end())
            {
                return false;
            }
        }
        return true;
    }

    bool hasCommonElement(const std::vector<SDP> &sdps)
    {   
        // Filter out the complementary sets
        std::vector<const SDP*> complementarySDPs;
        for (const auto &sdp : sdps)
        {
            if (sdp.isComplementary())
            {
                complementarySDPs.push_back(&sdp);
            }
        }

        // Check if there is a common element in the list of SDP sets
        for (size_t i = 0; i < complementarySDPs.size(); ++i)
        {
            for (size_t j = i + 1; j < complementarySDPs.size(); ++j)
            {
                // Find the intersection of the two sets
                std::vector<NodeID> commonElements;
                std::set_intersection(complementarySDPs[i]->begin(), complementarySDPs[i]->end(),
                                      complementarySDPs[j]->begin(), complementarySDPs[j]->end(),
                                      std::back_inserter(commonElements));

                // If there are common elements, return true
                if (!commonElements.empty())
                {
                    return true;
                }
            }
        }
        return false;
    }

    std::vector<Set> readPathsetsFromFile(const std::string &filename)
    {
        std::vector<Set> pathsets;
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filename);
        }
        std::string line;
        while (std::getline(file, line))
        {
            Set pathset;
            std::istringstream iss(line);
            NodeID node;
            while (iss >> node)
            {
                pathset.push_back(node);
            }
            if (!pathset.empty())
            {
                pathsets.push_back(std::move(pathset)); // Use move to avoid copying
            }
        }
        file.close();
        return pathsets;
    }
}