#pragma once
#include <pyrbd_plusplus/common.hpp>

namespace pyrbdpp::sdp
{
    // Declaration of the short types for the pathset module
    using PathSets = std::vector<Set>;
    using SDPSets = std::vector<SDP>;

    /**
     * @brief Eliminate the SPD set according to the algorithm from the paper "A simple algorithm for sum of disjoint products" by Ji Xing
     * Algorithm:
     * If the complementary sets have the same elements as non complementary sets, we can eliminate them from the complementary sets.
     * e.g. if we have the SPD set:
     *  A = {{2, 3, 1, 4, 6}}
     * -B = -{3, 4, 5}
     * 3 and 4 are in non-complementary set A, so we can eliminate them from the complementary set -B.
     * Result: {{2, 3, 1, 4, 6}, -{5}}
     *
     * @param sdpSets SPD sets with complementary and non-complementary sets
     * @return Eliminated SDP sets
     */
    SDPSets eliminateSDPSet(SDPSets &sdpSets);

    /**
     * @brief Absorbs the SDP set
     * If a SDP has a subset in the SDP set, it is absorbed by the SDP set, meaning remove the superset SDP.
     * e.g {3, 4}, {4, 5}, {3}, {5} -> {3}, {5}
     * 
     * @param sdpSets
     * @return absorbed SDP sets
     * @note The function will only absorb the complementary sets.
     */
    SDPSets absorbSDPSet(SDPSets sdpSets);

    /**
     * @brief Decomposes the SDP set according the algorithm from the paper "A simple algorithm for sum of disjoint products" by Ji Xing
     * Algorithm:
     * Assume the we have we SDP sets:
     *  A = {2, 7}
     * -B = -{3, 4, 5}
     * -C = -{1, 3, 7}
     * We found that in the complementary sets -B and -C we have the same elements, so we can decompose the SDP set into:
     * A (-{the same elements} + {the same elements} -{B without the same elements} -{C without the same elements})
     * e.g {2, 7} (-{3} + {3} -{4, 5} -{1, 7})
     * The result after expanding the add will be two SDPs:
     * SDP1: {2, 7}, -{3}
     * SDP2: {2, 7}, {3}, -{4, 5}, -{1, 7}
     *
     * @param sdpSets SDP set with complementary and non-complementary sets
     * @param results Vector of SDPSets to store the results of the decomposition
     * @return True if the rekursive decomposition was successful, false otherwise.
     * @note The function will recursively decompose the SDP sets until no more decompositions are possible.
     *       The result will be a vector of SDP sets, where each set is a decomposition of the original SDP set.
     *       The function will also eliminate and absorb the SDP sets before returning the result.
     *       This means that the result will not contain any redundant sets.
     */
    std::vector<SDPSets> decomposeSDPSet(SDPSets sdpSets);

    /**
     * @brief Sort the path sets acording to the algorithm from the paper "A Procedure for Generating the Sums of Disjoint Products" by Brijendra Singh
     * Part 1:
     * 1. Sort each set in the path sets in ascending order acording to the value e.g. {2, 4, 1} -> {1, 2, 4}, {3, 1, 7, 5} -> {1, 3, 5, 7}
     * 2. Sort the path sets in ascending order according to the size of the set e.g. {{1, 2, 4}, {1, 3, 5, 7}, {1, 2}} -> {{1, 2}, {1, 2, 4}, {1, 3, 5, 7}}
     * 3. Sort the sets with the same size acorrding to the increasing number of literals in commom with the preceding sets:
     *    e.g. for {{2, 7}, {1, 3, 7}, {1, 4, 6}, {2, 5, 6}, {1, 3, 5, 6}, {1, 4, 5, 7}, {2, 3, 4, 6}}
     *          1) {2, 7} is the first of the first group with size 2
     *          2) Preceding sets are: {2, 7}
     *             {1, 3, 7} has 1 literal in common with {2, 7}
     *             {1, 4, 6} has 0 literal in common with {2, 7}
     *             {2, 5, 6} has 1 literal in common with {2, 7}
     *             So the second group with size 3 is {{1, 4, 6}, {1, 3, 7}, {2, 5, 6}}
     *          3) Preceding sets are: {2, 7}, {1, 4, 6}, {1, 3, 7}, {2, 5, 6}
     *            {1, 3, 5, 6} has maximal 2 literals in common with {2, 5, 6}
     *            {1, 4, 5, 7} has maximal 2 literals in common with {1, 3, 7}
     *            {2, 3, 4, 6} has maximal 2 literals in common with {1, 4, 6}
     *             So the third group with size 4 is {{1, 3, 5, 6}, {1, 4, 5, 7}, {2, 3, 4, 6}}
     *          4) The final sorted pathSet is:
     *            {{2, 7}, {1, 4, 6}, {1, 3, 7}, {2, 5, 6}, {1, 3, 5, 6}, {1, 4, 5, 7}, {2, 3, 4, 6}}
     *
     * @param pathSets
     * @return Sorted pathSets
     * @note The pathSets are sorted in place, so the original pathSets will be modified.
     */
    PathSets sortPathSet(PathSets pathSets);

    /**
     * @brief Sums of Disjoint Products (SDP) algorithm from the papar "A Procedure for Gnerating the Sums of Disjoint Products" by Brijendra Singh
     * Algorithm:
     * Part 1:
     * See sortedPathSet() for the first part of the algorithm.
     *
     * Part 2:
     * 1. The final set is initalized with the first set in the sorted pathSet e.g. {{2, 7}}
     * 2. For each other set in the sorted pathSet make them disjoint with each previous sets:
     *    1) To make the sets disjoint, we add complementary elements to the existing sets in the final set,
     *      e.g. for the second set {1, 4, 6} it should be disjoint with first set {2, 7}.
     *      We define the RC list = {{2, 7} \ {1, 4, 6}} = {{2, 7}}, so we add  {2, 7}' to the second set.
     *      At the end we have the final set {{2, 7}, {2, 7}'{1, 4, 6}}.
     *      Be carefule {2, 7}' is not same as {2}'{7}', we evaluate the availability of {2, 7}' as 1 - ((avail of 2) * (avail of 7))
     *
     *    2) For the more complex case the RC list can be longer,
     *      e.g. for the third set {1, 3, 7} it should be disjoint with the first two sets {2, 7} and {1, 4, 6}.
     *      We define the RC list = {{{2, 7} \ {1, 3, 7}}, {1, 4, 6} \ {1, 3, 7}} = {{2}, {4, 6}}, so we add {2}'{4, 6}' to the third set.
     *      Now we have the final set {{2, 7}, {2, 7}'{1, 4, 6}, {2}'{4, 6}'{1, 3, 7}}.
     *      Be carefule the RC list can not have any super set of other set in the RC list,
     *      e.g. Example RC list = {{2}, {2, 7}, {3}, {3, 4, 5}} is not allowed, the super set must be removed so {{2}, {3}} is correct.
     *      This part is done by decomposeSDPSet() function.
     * 3. We add each set from pathset to the final set according to step 2.
     * 4. The final set is the result of the SDP algorithm.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param pathSets Path sets for the source and destination pair
     * @note The pathSets are sorted in place, so the original pathSets will be modified.
     *       The src and dst are not used in this function, but they are included for compatibility with other functions.  
     * @return Vector of SDP sets
     */
    std::vector<SDPSets> toSDPSet(NodeID src, NodeID dst, PathSets pathSets);

    /**
     * @brief The parallel version of the toSDPSet() function for large path sets.
     * The thread are used to parallelize the part2 of the SDP algorithm.
     * e.g if we have the path sets:
     * {{2, 7}, {1, 4, 6}, {1, 3, 7}, {2, 5, 6}, {1, 3, 5, 6}, {1, 4, 5, 7}, {2, 3, 4, 6}}
     * The function will create a thread for each set in the path sets and will process them in parallel.
     * Thread1: handle {1, 4, 6}
     * Thread2: handle {1, 3, 7}
     * Thread3: handle {2, 5, 6} usw.
     * The result will be saved unsorted in the final set, which will not affect the final availability.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param pathSets Path sets for the source and destination pair
     * @note The pathSets are sorted in place, so the original pathSets will be modified.
     *       The src and dst are not used in this function, but they are included for compatibility with other functions.  
     * @return Vector of SDP sets
     * @note This function is only used for the pathSets size larger equal than 1000 to avoid exceeding the number of available cores.
     *       If the path sets size smaller than 1000, the toSDPSet() function will be used.
     *       The function will be called automatically by the evalAvailParallel() function.
     */
    std::vector<SDPSets> toSDPSetParallel(NodeID src, NodeID dst, PathSets pathSets);

    /**
     * @brief Evaluate the availability of SDP set
     * @param probaMap ProbabilityMap containing the availability of each node
     * @param sdpSets SDP sets to evaluate the availability
     * @return Availability in double
     * @note The availability is directly calculated from the SDP sets.
     */
    double SDPSetToAvail(const ProbabilityMap &probaMap, const std::vector<SDPSets> &sdpSets);

    /**
     * @brief Evaluate the availability for a specific source and destination from a topology with SDP algorithm.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param probaMap ProbabilityMap containing the availability of each node
     * @param pathSets Path sets for the source and destination pair
     * @return Availability between source and destination in double
     */
    double evalAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, PathSets &pathSets);

    /**
     * @brief Evaluate the availability for a specific source and destination from a topology with SDP algorithm.
     * This function uses OpenMP to parallelize the evaluation of availability for the source and destination nodes.
     * @param src Source node ID
     * @param dst Destination node ID
     * @param probaMap ProbabilityMap containing the availability of each node
     * @param pathSets Path sets for the source and destination pair
     * @return Availability between source and destination in double
     * @note This function is only used for the pathSets size larger equal than 1000.
     * The function will not be used together with the evalAvailTopoParallel() function to avoid exceeding the number of available cores.
     * This function will call the toSDPSetParallel() function to convert the path sets to SDP sets.
     */
    double evalAvailParallel(NodeID src, NodeID dst, const ProbabilityMap &probaMap, PathSets &pathSets);

    /**
     * @brief Evaluate the availability for each pair of source and destination nodes in a topology with SDP algorithm.
     * @param nodePairs A vector of pairs of source and destination node IDs
     * @param probaMap Probability map containing the availability of each node
     * @param pathsetsList A vector of path sets for each node pair
     * @return List of (src, dst, availability) triples
     */
    std::vector<AvailTriple> evalAvailTopo(const NodePairs &nodePairs, const ProbabilityMap &probaMap, std::vector<PathSets> &pathsetsList);

    /**
     * @brief Evaluate the availability for each node pair in the nodePairs and topology file with SDP algorithm
     * This function uses OpenMP to parallelize the evaluation of availability for each pair of source and destination nodes.
     * @param nodePairs A vector of pairs of source and destination node IDs
     * @param probaMap Probability map containing the availability of each node
     * @param pathsetsList A vector of path sets for each node pair
     * @return List of (src, dst, availability) triples
     * @note This function can be used for every path sets size.
     * Do not use this function to call evalAvailParallel() function to keep the number of available cores not exceeded.
     */
    std::vector<AvailTriple> evalAvailTopoParallel(const NodePairs &nodePairs, const ProbabilityMap &probaMap, std::vector<PathSets> &pathsetsList);


} // namespace pyrbdpp::sdp