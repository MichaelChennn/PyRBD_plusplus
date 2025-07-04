#include <pyrbd_plusplus/sdp.hpp>
#include <pyrbd_plusplus/utils.hpp>
#include <numeric>
#include <chrono>
#include <mutex>

namespace pyrbdpp::sdp
{   

    using pyrbdpp::utils::isSubSet;

    SDPSets eliminateSDPSet(const SDPSets &sdpSets)
    {
        // Create a vector to store the elements that will be eliminated
        // from the non-complementary sets based on the complementary sets
        std::vector<int> eliminatedElements;

        for (const auto &SDP : sdpSets)
        {
            if (!SDP.isComplementary())
            {
                // If the set is non complementary, we add its elements to the eliminatedElements
                eliminatedElements.insert(eliminatedElements.end(), SDP.begin(), SDP.end());
            }
        }

        // Create a new set for the eliminated SDP sets
        SDPSets eliminatedSet;

        for (const auto &SDP : sdpSets)
        {
            if (SDP.isComplementary())
            {
                // If the set is complementary, check if it has any elements in common with the eliminatedElements
                std::vector<int> newSet;
                std::set_difference(SDP.begin(), SDP.end(), eliminatedElements.begin(), eliminatedElements.end(),
                                    std::back_inserter(newSet));
                eliminatedSet.emplace_back(true, newSet);
            }
            else
            {
                // If the set is non complementary, we just add it to the eliminatedSet
                eliminatedSet.push_back(SDP);
            }
        }

        return eliminatedSet;
    }

    SDPSets absorbSDPSet(const SDPSets &sdpSets)
    {   
        // Initialize a vector to store the absorbed SDPs
        SDPSets absorbedSDPs;

        // If the set is complementary, add it to the complementarySDPs else add it to absorbedSDPs
        SDPSets complementarySDPs;
        SDPSets noncomplementarySDPs;

        for (const auto &SDP : sdpSets)
        {
            if (SDP.isComplementary())
            {
                complementarySDPs.push_back(SDP);
            }
            else
            {
                absorbedSDPs.push_back(SDP);
            }
        }

        std::vector<bool> absorbed(complementarySDPs.size(), false);

        for (size_t i = 0; i < complementarySDPs.size(); ++i)
        {
            if (absorbed[i])
                continue; // Skip if already absorbed

            const auto &currentSDP = complementarySDPs[i];

            for (size_t j = 0; j < complementarySDPs.size(); ++j)
            {
                if (i == j || absorbed[j])
                    continue; // Skip self and already absorbed sets

                const auto &otherSDP = complementarySDPs[j];

                // Check if currentSDP is a subset of otherSDP
                if (isSubSet(currentSDP.getSet(), otherSDP.getSet()))
                {
                    absorbed[j] = true; // Mark the superset as absorbed
                }

                if (isSubSet(otherSDP.getSet(), currentSDP.getSet()))
                {
                    absorbed[i] = true; // Mark the subset as absorbed
                }
            }

            if (!absorbed[i])
            {
                absorbedSDPs.push_back(currentSDP); // Keep the current set if not absorbed
            }
        }

        return absorbedSDPs;
    }

    bool decomposeSDPSet(const SDPSets &sdpSets, std::vector<SDPSets> &results)
    {
        // Find the normal sets and the complementary sets
        SDPSets complementarySDPs;

        for (const auto &SDP : sdpSets)
        {
            if (SDP.isComplementary())
            {
                // If the set is complementary, add it to the complementarySDPs
                complementarySDPs.push_back(SDP);
            }
        }

        // If there are no complementary sets, return the original SDPs
        if (complementarySDPs.empty())
        {
            results.push_back(sdpSets);
            return true;
        }

        std::pair<SDP, SDP> commonSDPpair;
        std::vector<int> commonElements;

        bool firstCommonPairFound = false;
        // Find only the first common elements in the complementary sets and handle it
        for (size_t i = 0; i < complementarySDPs.size(); ++i)
        {
            const auto &compSDP = complementarySDPs[i];

            // Check if the set has any elements in common with other complementary sets
            for (size_t j = i + 1; j < complementarySDPs.size(); ++j)
            {
                const auto &otherCompSDP = complementarySDPs[j];

                // Find the intersection of the two sets
                std::set_intersection(compSDP.begin(), compSDP.end(),
                                      otherCompSDP.begin(), otherCompSDP.end(),
                                      std::back_inserter(commonElements));

                // If there are common elements, add them to the commonElements
                if (!commonElements.empty())
                {
                    commonSDPpair = std::make_pair(compSDP, otherCompSDP);
                    firstCommonPairFound = true;
                    break;
                }
            }

            // If we found the first common pair, break the loop
            if (firstCommonPairFound)
            {
                break;
            }
        }

        // If there is no more common elements, add the original SDP set to the results and return true
        if (commonElements.empty())
        {
            // Add the original SDPs to the results
            results.push_back(sdpSets);
            return true;
        }

        // Add all SDP except common pair to the normal SDP set
        SDPSets normalSDPs;
        for (const auto &SDP : sdpSets)
        {
            if (!SDP.equals(commonSDPpair.first) && !SDP.equals(commonSDPpair.second))
            {
                normalSDPs.push_back(SDP);
            }
        }

        // Remove the common elements from the common pairs
        for (const auto &elem : commonElements)
        {
            commonSDPpair.first.remove(elem);
            commonSDPpair.second.remove(elem);
        }

        // Create two new SDP sets for the decomposed sets
        SDPSets decomposedSDPs1;
        SDPSets decomposedSDPs2;

        // For the first new SDP sets we add all normal SDPs and the complementary SDP with the common elements
        SDP commonSDP1(true, commonElements);
        decomposedSDPs1.insert(decomposedSDPs1.end(), normalSDPs.begin(), normalSDPs.end());
        decomposedSDPs1.push_back(commonSDP1);

        // For the second new SDP sets we add
        // 1. all normal SDP
        // 2. Non complementary SDP with the common elements
        // 3. Complementary SDP with elements from the first common pair without common elements
        // 4. Complementary SDP with elements from the second common pair without common elements
        SDP commonSDP2(false, commonElements);
        decomposedSDPs2.insert(decomposedSDPs2.end(), normalSDPs.begin(), normalSDPs.end());
        decomposedSDPs2.push_back(commonSDP2);
        decomposedSDPs2.push_back(commonSDPpair.first);
        decomposedSDPs2.push_back(commonSDPpair.second);

        // Eliminate the SDP sets to remove any redundant sets
        decomposedSDPs1 = eliminateSDPSet(decomposedSDPs1);
        decomposedSDPs2 = eliminateSDPSet(decomposedSDPs2);

        // Absorb the SDP sets to remove any redundant sets
        decomposedSDPs1 = absorbSDPSet(decomposedSDPs1);
        decomposedSDPs2 = absorbSDPSet(decomposedSDPs2);

        // recursively decompose new SDPs in results
        return decomposeSDPSet(decomposedSDPs1, results) && decomposeSDPSet(decomposedSDPs2, results);
    }

    PathSets sortPathSet(PathSets &pathSets)
    {
        // Sort the number in each set of pathSets in ascending order
        for (auto &set : pathSets)
        {
            std::sort(set.begin(), set.end());
        }

        // Sort the pathSets in ascending order according to the size of the set and the number in each set
        std::sort(pathSets.begin(), pathSets.end(), [](const std::vector<int> &a, const std::vector<int> &b)
                  {
            if (a.size() != b.size())
            {
                return a.size() < b.size();
            }
            return a < b; });

        // Separate the pathSets into a map according to the size of the set
        std::map<int, PathSets> pathSetMap;
        for (const auto &set : pathSets)
        {
            pathSetMap[set.size()].push_back(set);
        }

        // Save the sorted pathSet
        std::vector<std::vector<int>> sortedPathSet;

        // Reserve space for the sortedPathSet to avoid reallocation
        sortedPathSet.reserve(pathSets.size());

        // Sort the sets with the same size according to the increasing maximal number of literals in common with the preceding sets
        for (auto &pair : pathSetMap)
        {
            // pair.first is the size of the set, pair.second is the vector of sets with this size
            auto &unsortedSets = pair.second;

            // If the sortedPathSet is empty, just add all sets with this size
            if (sortedPathSet.empty())
            {
                sortedPathSet.insert(sortedPathSet.end(), unsortedSets.begin(), unsortedSets.end());
                continue;
            }

            // Create a vector to store the maximal number of literals in common with the preceding sets
            Set maxCommonLiteralsLst;

            // For each set stores the maximal number of literals in common with the preceding sets
            maxCommonLiteralsLst.reserve(unsortedSets.size());
            std::transform(unsortedSets.begin(), unsortedSets.end(), std::back_inserter(maxCommonLiteralsLst), [&sortedPathSet](const std::vector<int> &set)
                           {
                int maxCommonLiterals = 0;
                for (const auto &precedSet : sortedPathSet)
                {
                    // Count the number of literals in common
                    int commonLiterals = std::count_if(set.begin(), set.end(), [&precedSet](int elem) {
                        return std::find(precedSet.begin(), precedSet.end(), elem) != precedSet.end();
                    });

                    // If the current set has more literals in common, update the best set
                    if (commonLiterals > maxCommonLiterals)
                    {
                        maxCommonLiterals = commonLiterals;
                    }
                }
                return maxCommonLiterals; });

            // Create new indices for the unsortedSets based on the maximal number of literals in common in ascending order
            std::vector<size_t> indices(unsortedSets.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&maxCommonLiteralsLst](size_t a, size_t b)
                      { return maxCommonLiteralsLst[a] < maxCommonLiteralsLst[b]; });

            // Add the sorted sets to the sortedPathSet based on the indices
            for (size_t idx : indices)
            {
                sortedPathSet.push_back(unsortedSets[idx]);
            }
        }

        // Now we have the sorted pathSet, we can return it
        return sortedPathSet;
    }

    std::vector<SDPSets> toSDPSet(NodeID src, NodeID dst, PathSets &pathSets)
    {
        // Sort the pathSets
        std::vector<std::vector<int>> sortedPathSet = sortPathSet(pathSets);

        // Now we have the sorted pathSet, we can apply the second part of the SDP algorithm

        // Use the SDP set to store the final result, initialized with the first set in the sorted pathSet
        SDPSets firstSDP;
        firstSDP.emplace_back(false, sortedPathSet.front());
        std::vector<SDPSets> finalSDP = {{firstSDP}};

        // Iterate over the sorted pathSet starting from the second set and find only the first pair of sets that have common elements
        for (size_t i = 1; i < sortedPathSet.size(); ++i)
        {
            // Get the current set
            const auto &currentSet = sortedPathSet[i];

            // Initialize the final SDP set for the current set
            SDPSets currentSDPs;
            currentSDPs.emplace_back(false, currentSet);

            // Iterate over the previous sets in sortedPathSet
            for (size_t j = 0; j < i; ++j)
            {
                const auto &precedingSet = sortedPathSet[j];

                // Create the RC set: elements in precedingSet but not in currentSet
                std::vector<int> RC;
                std::set_difference(precedingSet.begin(), precedingSet.end(),
                                    currentSet.begin(), currentSet.end(),
                                    std::back_inserter(RC));

                // Add not empty RC sets to the currentSDPs
                if (!RC.empty())
                {
                    currentSDPs.emplace_back(true, RC);
                }
            }

            // Absorb the currentSDPs to remove any redundant sets
            currentSDPs = absorbSDPSet(currentSDPs);

            // Decompose the currentSDPs and add all results to the finalSDP
            std::vector<SDPSets> decomposedResults;
            if (decomposeSDPSet(currentSDPs, decomposedResults))
            {
                finalSDP.insert(finalSDP.end(), decomposedResults.begin(), decomposedResults.end());
            }
            else
            {
                // If decomposition fails, Error
                std::cerr << "Error: Decomposition failed for currentSDPs." << std::endl;
            }
        }

        return finalSDP;
    }

    std::vector<SDPSets> toSDPSetParallel(NodeID src, NodeID dst, PathSets &pathSets)
    {   
        // Check if the pathSets size is less than 1000
        if (pathSets.size() < 1000)
        {
            // If the pathSets size is less than 1000, use the single-threaded version
            std::cerr << "Warning: Trying to use multithreaded version of toSDPSet for small pathSets size." << std::endl;
            std::cerr << "Using single-threaded version instead." << std::endl;
            // Call the single-threaded version of toSDPSet
            return toSDPSet(src, dst, pathSets);
        }
        
        // Sort the pathSets
        std::vector<std::vector<int>> sortedPathSet = sortPathSet(pathSets);

        // Now we have the sorted pathSet, we can apply the second part of the SDP algorithm

        // Initialize a the final thread results
        std::vector<std::vector<SDPSets>> threadResults(sortedPathSet.size());

        // Use the SDP set to store the final result, initialized with the first set in the sorted pathSet
        SDPSets firstSDP;
        firstSDP.emplace_back(false, sortedPathSet.front());
        threadResults[0] = {firstSDP};

        // Iterate over the sorted pathSet starting from the second set and find only the first pair of sets that have common elements
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 1; i < sortedPathSet.size(); ++i)
        {
            // Get the current set
            const auto &currentSet = sortedPathSet[i];

            // Initialize the final SDP set for the current set
            SDPSets currentSDPs;
            currentSDPs.emplace_back(false, currentSet);

            // Iterate over the previous sets in sortedPathSet
            for (size_t j = 0; j < i; ++j)
            {
                const auto &precedingSet = sortedPathSet[j];

                // Create the RC set: elements in precedingSet but not in currentSet
                std::vector<int> RC;
                std::set_difference(precedingSet.begin(), precedingSet.end(),
                                    currentSet.begin(), currentSet.end(),
                                    std::back_inserter(RC));

                // Add not empty RC sets to the currentSDPs
                if (!RC.empty())
                {
                    currentSDPs.emplace_back(true, RC);
                }
            }

            // Absorb the currentSDPs to remove any redundant sets
            currentSDPs = absorbSDPSet(currentSDPs);

            // Decompose the currentSDPs and add all results to the threadResults
            if (decomposeSDPSet(currentSDPs, threadResults[i])){ }
            else
            {
                // If decomposition fails, Error
                std::cerr << "Error: Decomposition failed for currentSDPs." << std::endl;
            }
        }

        // Add the thread results to the finalSDP
        std::vector<SDPSets> finalSDP;
        for (const auto &result : threadResults)
        {
            finalSDP.insert(finalSDP.end(), result.begin(), result.end());
        }

        return finalSDP;
    }

    DebugInfo toSDPSetDebug(NodeID src, NodeID dst, PathSets &pathSets)
    {   
        // Initialize DebugInfo to store the debug information
        DebugInfo debugInfo;

        // Sort the pathSets
        std::vector<std::vector<int>> sortedPathSet = sortPathSet(pathSets);

        // Now we have the sorted pathSet, we can apply the second part of the SDP algorithm

        // Use the SDP set to store the final result, initialized with the first set in the sorted pathSet
        SDPSets firstSDP;
        firstSDP.emplace_back(false, sortedPathSet.front());
        std::vector<SDPSets> finalSDP = {{firstSDP}};

        // Iterate over the sorted pathSet starting from the second set and find only the first pair of sets that have common elements
        for (size_t i = 1; i < sortedPathSet.size(); ++i)
        {   
            // Measure the time taken for each iteration
            auto start = std::chrono::high_resolution_clock::now();

            // Get the current set
            const auto &currentSet = sortedPathSet[i];

            // Initialize the final SDP set for the current set
            SDPSets currentSDPs;
            currentSDPs.emplace_back(false, currentSet);

            // Iterate over the previous sets in sortedPathSet
            for (size_t j = 0; j < i; ++j)
            {
                const auto &precedingSet = sortedPathSet[j];

                // Create the RC set: elements in precedingSet but not in currentSet
                std::vector<int> RC;
                std::set_difference(precedingSet.begin(), precedingSet.end(),
                                    currentSet.begin(), currentSet.end(),
                                    std::back_inserter(RC));

                // Add not empty RC sets to the currentSDPs
                if (!RC.empty())
                {
                    currentSDPs.emplace_back(true, RC);
                }
            }

            // Absorb the currentSDPs to remove any redundant sets
            currentSDPs = absorbSDPSet(currentSDPs);

            // Decompose the currentSDPs and add all results to the finalSDP
            std::vector<SDPSets> decomposedResults;
            if (decomposeSDPSet(currentSDPs, decomposedResults))
            {
                finalSDP.insert(finalSDP.end(), decomposedResults.begin(), decomposedResults.end());
            }
            else
            {
                // If decomposition fails, Error
                std::cerr << "Error: Decomposition failed for currentSDPs." << std::endl;
            }

            // Measure the time taken for the current iteration
            auto end = std::chrono::high_resolution_clock::now();

            // Calculate the duration in milliseconds
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            debugInfo[i] = {finalSDP.size(), duration.count()};
        }

        return debugInfo;
    }

    double SDPSetToAvail(const ProbabilityMap &probaMap, const std::vector<SDPSets> &sdpSets)
    {
        double availability = 0.0;

        // Iterate over each set in the SDP set
        for (const auto &set : sdpSets)
        {
            double setAvailability = 1.0;

            // Compute the availability for the current set
            for (const auto &SDP : set)
            {

                if (SDP.isComplementary())
                {
                    double tmp_avail = 1.0;

                    // If the set is complementary, use the negative probability
                    for (const auto &elem : SDP.getSet())
                    {
                        tmp_avail *= probaMap[elem];
                    }

                    setAvailability *= (1.0 - tmp_avail);
                }
                else
                {
                    // If the set is normal, use the positive probability
                    for (const auto &elem : SDP.getSet())
                    {
                        setAvailability *= probaMap[elem];
                    }
                }
            }

            // Add the set availability to the total availability
            availability += setAvailability;
        }

        return availability;
    }

    double evalAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, PathSets &pathSets)
    {
        // Convert the pathSets to the SDP set
        std::vector<SDPSets> SDPs = toSDPSet(src, dst, pathSets);

        // Calculate the availability of the SDP set
        double availability = SDPSetToAvail(probaMap, SDPs);

        return availability;
    }

    double evalAvailParallel(NodeID src, NodeID dst, const ProbabilityMap &probaMap, PathSets &pathSets)
    {
        // Convert the pathSets to the SDP set
        std::vector<SDPSets> SDPs = toSDPSetParallel(src, dst, pathSets);

        // Calculate the availability of the SDP set
        double availability = SDPSetToAvail(probaMap, SDPs);

        return availability;
    }

    std::vector<AvailTriple> evalAvailTopo(const NodePairs &nodePairs, const ProbabilityMap &probaMap, std::vector<PathSets> &pathsetsList)
    {
        std::vector<AvailTriple> availList;

        for (size_t i = 0; i < nodePairs.size(); ++i)
        {   
            const auto &[src, dst] = nodePairs[i];
            auto &pathSets = pathsetsList[i];

            double availability = evalAvail(src, dst, probaMap, pathSets);
            availList.emplace_back(src, dst, availability);
        }

        return availList;
    }

    std::vector<AvailTriple> evalAvailTopoParallel(const NodePairs &nodePairs, const ProbabilityMap &probaMap, std::vector<PathSets> &pathsetsList)
    {
        std::vector<AvailTriple> availList(nodePairs.size());

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < nodePairs.size(); ++i)
        {
            const auto &[src, dst] = nodePairs[i];
            auto &pathSets = pathsetsList[i];

            double availability = evalAvail(src, dst, probaMap, pathSets);
            availList[i] = std::make_tuple(src, dst, availability);
        }

        return availList;
    }

} // namespace pyrbdpp::sdp
