#include <pyrbd_plusplus/mcs.hpp>
#include <algorithm>
#include <chrono>
#include <omp.h>

namespace pyrbdpp::mcs
{   

    ProbaSets toProbaSet(NodeID src, NodeID dst, MinCutSets minCutSets)
    {

        // remove the sets: {src} and {dst}
        minCutSets.erase(std::remove(minCutSets.begin(), minCutSets.end(), std::vector<int>{src}), minCutSets.end());
        minCutSets.erase(std::remove(minCutSets.begin(), minCutSets.end(), std::vector<int>{dst}), minCutSets.end());

        // Check if the minCutSets is empty after removing the src and dst sets
        if (minCutSets.empty())
            return {};

        // Inverse the minimal cut sets
        for (auto &set : minCutSets)
        {
            std::transform(set.begin(), set.end(), set.begin(), [](int x)
                           { return -x; });
        }

        // Initialize the probability sets
        ProbaSets probaSets;

        // Reserve space for the probability sets
        probaSets.reserve(minCutSets.size() * 3); // Assuming each min cut set can generate up to 3 disjoint sets

        // convert the minimal cut sets to the probability sets
        while (minCutSets.size() > 0)
        {
            if (minCutSets.size() == 1)
            {
                probaSets.push_back(minCutSets[0]);
                break;
            }

            Set selectedSet = minCutSets.front();

            probaSets.push_back(selectedSet);

            std::vector<std::vector<int>> remainingSets(minCutSets.begin() + 1, minCutSets.end());

            minCutSets.clear();

            for (const auto &set : remainingSets)
            {
                // Create the disjoint sets from the selectedSet and each set in the remaining sets
                DisjointSets disjointSets = makeDisjointSet(selectedSet, set);

                // Insert the disjoint sets to the minCutSet
                minCutSets.insert(minCutSets.end(), disjointSets.begin(), disjointSets.end());
            }
        }

        return probaSets;
    }

    DebugInfo toProbaSetDebug(NodeID src, NodeID dst, MinCutSets minCutSets)
    {
        // Initialize the debug information
        DebugInfo debugInfo;

        // remove the sets: {src} and {dst}
        minCutSets.erase(std::remove(minCutSets.begin(), minCutSets.end(), std::vector<int>{src}), minCutSets.end());
        minCutSets.erase(std::remove(minCutSets.begin(), minCutSets.end(), std::vector<int>{dst}), minCutSets.end());

        // Check if the minCutSets is empty after removing the src and dst sets
        if (minCutSets.empty())
            return debugInfo;

        // Inverse the minimal cut sets
        for (auto &set : minCutSets)
        {
            std::transform(set.begin(), set.end(), set.begin(), [](int x)
                           { return -x; });
        }

        // Initialize the probability sets
        ProbaSets probaSets;

        // Reserve space for the probability sets
        probaSets.reserve(minCutSets.size() * 3); // Assuming each min cut set can generate up to 3 disjoint sets

        // Initialize the iteration counter
        int iteration = 0;

        // convert the minimal cut sets to the probability sets
        while (minCutSets.size() > 0)
        {
            // Measure the time taken for each iteration
            auto start = std::chrono::high_resolution_clock::now();

            if (minCutSets.size() == 1)
            {
                probaSets.push_back(minCutSets[0]);
                break;
            }

            const Set &selectedSet = minCutSets.front();

            probaSets.push_back(selectedSet);

            std::vector<std::vector<int>> remainingSets(minCutSets.begin() + 1, minCutSets.end());

            minCutSets.clear();

            for (const auto &set : remainingSets)
            {
                // Create the disjoint sets from the selectedSet and each set in the remaining sets
                DisjointSets disjoint_sets = makeDisjointSet(selectedSet, set);

                // Insert the disjoint sets to the minCutSet
                minCutSets.insert(minCutSets.end(), disjoint_sets.begin(), disjoint_sets.end());
            }

            // Measure the end time
            auto end = std::chrono::high_resolution_clock::now();

            // Calculate the duration of the iteration
            std::chrono::duration<double> duration = end - start;

            // Save the debug information
            debugInfo[iteration] = {probaSets.size(), duration.count()};

            // Increment the iteration counter
            iteration++;
        }

        return debugInfo;
    }

    double probaSetToAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const ProbaSets &probaSet)
    {
        // Save the final result
        double unavil = 0.0;
        for (const auto &set : probaSet)
        {
            // Save the temp result
            double temp = 1.0;
            for (const auto &num : set)
            {
                temp *= probaMap[num];
            }
            unavil += temp;
        }

        double avail = 1.0 - unavil;
        double result = probaMap[src] * probaMap[dst] * avail;

        return result;
    }

    
    double evalAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const MinCutSets &minCutSets)
    {
        // Convert the minimal cut sets to the probability sets
        ProbaSets probaSets = toProbaSet(src, dst, minCutSets);

        // Compute the availability
        return probaSetToAvail(src, dst, probaMap, probaSets);
    }

    std::vector<AvailTriple> evalAvailTopo(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<MinCutSets> &minCutSetsList)
    {
        std::vector<AvailTriple> availList;

        for (size_t i = 0; i < nodePairs.size(); ++i)
        {   
            const auto &[src, dst] = nodePairs[i];
            const auto &minCutSets = minCutSetsList[i];

            double availability = evalAvail(src, dst, probaMap, minCutSets);
            availList.emplace_back(src, dst, availability);
        }

        return availList;
    }

    std::vector<AvailTriple> evalAvailTopoParallel(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<MinCutSets> &minCutSetsList)
    {
        std::vector<AvailTriple> availList(nodePairs.size());

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < nodePairs.size(); ++i)
        {
            const auto &[src, dst] = nodePairs[i];
            const auto &minCutSets = minCutSetsList[i];

            double availability = evalAvail(src, dst, probaMap, minCutSets);
            availList[i] = std::make_tuple(src, dst, availability);
        }

        return availList;
    }

} // namespace pyrbdpp::mcs
