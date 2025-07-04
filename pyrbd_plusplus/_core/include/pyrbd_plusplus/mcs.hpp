#pragma once
#include <pyrbd_plusplus/common.hpp>


namespace pyrbdpp::mcs
{
    // Declaration of the short types for the MCS module
    using MinCutSets   = std::vector<Set>;
   
    /**
     * @brief Convert the minimal cut sets to the probability sets
     * Algorithm:
     * 1. Remove the set: {src} and {dst} from the minimal cut sets
     * 2. Inverse the minimal cut sets, e.g. minCutSets = {{1, 2, 3}, {1, 2, -3}} -> {{-1, -2, -3}, {-1, -2, 3}}
     * 3. Select the left most set as the selected_set and the remaining minimal cut sets as the right remaining sets 
     *    e.g minimal cut sets = {set1, set2, set3}, selected_set = set1, remaining sets = {set2, set3}
     * 4. Create the disjoint sets from the set1 and each set in the remaining sets e.g. makeDisjointSet(set1, set2), makeDisjointSet(set1, set3)
     * 5. The new disjoint set should follow the rules:
     *    Absorption: x + xy = x
     * 6. Add the new disjoint sets to the probability sets
     * 7. repeat the steps 2-5 for the second left most set and the remaining sets are the right remaining sets until all set are selected once
     *    e.g. selected_set = set2, remaining sets = {set3}
     * @param src Source node ID
     * @param dst Destination node ID
     * @param minCutSets Minimal cut sets for the source and destination pair
     * @return Probability sets
     */
    ProbaSets toProbaSet(NodeID src, NodeID dst, MinCutSets minCutSets);


    /**
     * @brief Debug version of toProbaSet
     * This function is used to debug the toProbaSet function.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param minCutSets Minimal cut sets for the source and destination pair
     * @note The debug information is saved in a map with the iteration as key and the current ProbaSet size and duration as value.
     *       The duration is the time taken to execute one iteration of the algorithm.
     *       The map is used to track the progress of the algorithm and to debug the performance
     *       of the algorithm.
     * @return Debug information: {iteration : (current ProbaSet size, duration)
     */
    DebugInfo toProbaSetDebug(NodeID src, NodeID dst, MinCutSets minCutSets);

    
    /**
     * @brief Compute the probability of the given probability set and the probability map
     * e.g. probability sets = {{1, 2, 3}, {1, 2, -3}}, probability map = {1:0.1, 2:0.2, 3:0.3}
     * then the probability is 0.1 * 0.2 * 0.3 + 0.1 * 0.2 * 0.7 = 0.02
     * @param src Source node ID
     * @param dst Destination node ID
     * @param probabilityMap Probability map containing the availability probabilities for each node
     * @param probabilitySets Probability sets containing the sets of nodes
     * @return Availability between source and destination in double
     * @note The result from toProbaSet is used to calculate the unavailability.
     *       So the final result calculates (1 - unavailability) multiplies the source and destination nodes' probabilities.
     */
    double probaSetToAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const ProbaSets &probaSets);

    /**
     * @brief Evaluate the availability for a specific source and destination from a topology.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param probMap Probability map containing the availability probabilities for each node
     * @param minCutSets Minimal cut sets for the source and destination pair
     * @return Availability between source and destination in double
     */
    double evalAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const MinCutSets &minCutSets);

    /**
     * @brief Evaluate the availability for each pair of source and destination nodes in a topology.
     * @param nodePairs A vector of pairs of source and destination node IDs
     * @param probaMap Probability map containing the availability probabilities for each node
     * @param minCutSetsList A vector of minimal cut sets for each pair of source and destination nodes
     * @return List of (src, dst, availability) tuples
     */
    std::vector<AvailTriple> evalAvailTopo(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<MinCutSets> &minCutSetsList);

    /**
     * @brief Evaluate the availability for each pair of source and destination nodes in a topology in parallel.
     * This function uses OpenMP to parallelize the evaluation of availability for each pair of source and destination nodes.
     * Threads are created to compute the availability for each pair concurrently, improving performance for large topologies.
     * @param nodePairs A vector of pairs of source and destination node IDs
     * @param probaMap Probability map containing the availability probabilities for each node
     * @param minCutSetsList A vector of minimal cut sets for each pair of source and destination nodes
     * @return List of (src, dst, availability) tuples
     */
    std::vector<AvailTriple> evalAvailTopoParallel(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<MinCutSets> &minCutSetsList);


} // namespace pyrbdpp::mcs