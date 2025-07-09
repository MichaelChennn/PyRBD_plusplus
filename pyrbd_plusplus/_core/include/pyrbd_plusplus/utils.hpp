#pragma once
#include <pyrbd_plusplus/common.hpp>
#include <algorithm>

namespace pyrbdpp::utils
{   
    /**
     * @brief Check if set1 is a subset of set2.
     * @param set1 The potential subset.
     * @param set2 The set to check against.
     * @return True if set1 is a subset of set2, false otherwise.
     */
    bool isSubSet(const Set &set1, const Set &set2);

    /**
     * @brief Check if there is a common element in a list of complementary SDP sets
     * @param SDPs The list of SDP sets to check.
     * @return True if there is a common element, false otherwise.
     */
    bool hasCommonElement(const std::vector<SDP> &sdps);

    /**
     * @brief Read pathsets from a file and generate a vector of sets.
     * @param filename The name of the file containing pathsets.
     * @return A vector of sets representing the pathsets.
     * @note Each line in the file should contain space-separated integers representing a pathset.
     */
    std::vector<Set> readPathsetsFromFile(const std::string &filename);
}