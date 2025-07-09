#pragma once
#include <pyrbd_plusplus/common.hpp>
#include <pyrbd_plusplus/sdp.hpp>
#include <algorithm>

namespace pyrbdpp::utils
{   
    // Declaration of the short types for the pathset module
    using PathSets = pyrbdpp::sdp::PathSets; // Only visible in pyrbdpp::sdp namespace
    using SDPSets = pyrbdpp::sdp::SDPSets; // Only visible in pyrbdpp::sdp namespace

    /**
     * @brief Check if sdp1 is a subset of sdp2.
     * @param sdp1 The potential subset.
     * @param sdp2 The set to check against.
     * @return True if sdp1 is a subset of sdp2, false otherwise.
     */
    bool isSubSet(const SDP &sdp1, const SDP &sdp2);

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

    /**
     * @brief Write SDPSets to file
     * @param sdpSets The SDPSets to write
     * @param filename The output filename (default: "SDPSets.txt")
     */
    void writeSDPSetsToFile(const std::vector<std::vector<SDP>> &sdpSets, const std::string &filename = "SDPSets.txt");

    /**
     * @brief Convert a Set to string representation
     * @param set The set to convert
     * @return String representation of the set
     */
    std::string toString(const Set &set);

     /**
     * @brief Convert SDPSets to string representation
     * @param sdpSets The SDPSets to convert
     * @return String representation of the SDPSets
     */
    std::string toString(const SDPSets &sdpSets);

    /**
     * @brief Convert vector of SDPSets to string representation
     * @param vectorSDPSets The vector of SDPSets to convert
     * @return String representation of the vector of SDPSets
     */
    std::string toString(const std::vector<SDPSets> &vectorSDPSets);

     /**
     * @brief Convert PathSets to string representation
     * @param pathSets The pathSets to convert
     * @return String representation of the pathSets
     */
    std::string toString(const PathSets &pathSets);
}