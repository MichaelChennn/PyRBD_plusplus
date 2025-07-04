#pragma once
#include <pyrbd_plusplus/common.hpp>

namespace pyrbdpp::pathset
{
    // Declaration of the short types for the pathset module
    using PathSets = std::vector<Set>;

    /**
     * @brief Convert the path sets to the probability sets
     * Algorithm:
     * 1. Select the left most set as the selected_set and the remaining path sets as the right remaining sets e.g path sets = {set1, set2, set3}
     * then selected_set = set1, remaining sets = {set2, set3}
     * 2. Create the disjoint sets from the set1 and each set in the remaining sets e.g. makeDisjointSet(set1, set2), makeDisjointSet(set1, set3)
     * 4. The new disjoint set should follow the rules:
     * Absorption: x + xy = x
     * 5. Add the new disjoint sets to the probability sets
     * 6. repeat the steps 2-5 for the second left most set and the remaining sets are the right remaining sets e.g. selected_set = set2, remaining sets = {set3}
     * @param src Source node ID
     * @param dst Destination node ID
     * @param pathSets Path sets for the source and destination pair
     * @return Probability sets
     */
    ProbaSets toProbaSet(NodeID src, NodeID dst, PathSets pathSets);

    /**
     * @brief Debug version of toProbaSet
     * This function is used to debug the toProbaSet function.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param pathSets Path sets for the source and destination pair
     * @note The debug information is saved in a map with the iteration as key and the current ProbaSet size and duration as value.
     *       The duration is the time taken to execute one iteration of the algorithm.
     *       The map is used to track the progress of the algorithm and to debug the performance
     *       of the algorithm.
     * @return Debug information: {iteration : (current ProbaSet size, duration)}
     */
    DebugInfo toProbaSetDebug(NodeID src, NodeID dst, PathSets pathSets);

    
    /**
     * @brief Compute the probability of the given probability set and the probability map
     * e.g. probability sets = {{1, 2, 3}, {1, 2, -3}}, probability map = {1:0.1, 2:0.2, 3:0.3}
     * then the probability is 0.1 * 0.2 * 0.3 + 0.1 * 0.2 * 0.7 = 0.02
     * @param src Source node ID
     * @param dst Destination node ID
     * @param probabilityMap Probability map containing the availability probabilities for each node
     * @param probabilitySets Probability sets containing the sets of nodes
     * @return Availability between source and destination in double
     * @note The availability do not multiply the source and destination nodes' probabilities at the end.
     *       So the final result is directly calculated from the probability sets.
     */
    double probaSetToAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const ProbaSets &probaSets);

    /**
     * @brief Evaluate the availability for a specific source and destination from a topology.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param probMap Probability map containing the availability probabilities for each node
     * @param pathSets Path sets for the source and destination pair
     * @return Availability between source and destination in double
     */
    double evalAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const PathSets &pathSets);

    /**
     * @brief Evaluate the availability for each pair of source and destination nodes in a topology.
     * @param nodePairs A vector of pairs of source and destination node IDs
     * @param probaMap Probability map containing the availability probabilities for each node
     * @param pathsetsList A vector of path sets for each pair of source and destination nodes
     * @return List of (src, dst, availability) tuples
     */
    std::vector<AvailTriple> evalAvailTopo(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<PathSets> &pathsetsList);

    /**
     * @brief Evaluate the availability for each pair of source and destination nodes in a topology in parallel.
     * This function uses OpenMP to parallelize the evaluation of availability for each pair of source and destination nodes.
     * @param nodePairs A vector of pairs of source and destination node IDs
     * @param probaMap Probability map containing the availability probabilities for each node
     * @param pathsetsList A vector of path sets for each pair of source and destination nodes
     * @return List of (src, dst, availability) tuples
     */
    std::vector<AvailTriple> evalAvailTopoParallel(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<PathSets> &pathsetsList);

} // namespace pyrbdpp::pathset