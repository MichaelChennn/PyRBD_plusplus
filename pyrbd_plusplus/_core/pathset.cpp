#include <pyrbd_plusplus/pathset.hpp>
#include <algorithm>
#include <chrono>
#include <omp.h>

namespace pyrbdpp::pathset
{

    ProbaSets toProbaSet(NodeID src, NodeID dst, PathSets pathSets)
    {   
        // Check if the pathSets is empty
        if (pathSets.empty())
        {
            return {};
        }
        
        // Initialize the probability sets
        ProbaSets probaSets;

        // Reserve space for the probability sets
        probaSets.reserve(pathSets.size() * 3); // Assuming each path set can generate up to 3 disjoint sets

        // convert the path sets to the probability sets
        while (pathSets.size() > 0)
        {
            if (pathSets.size() == 1)
            {
                probaSets.push_back(pathSets[0]);
                break;
            }

            Set selectedSet = pathSets.front();

            probaSets.push_back(selectedSet);

            std::vector<std::vector<int>> remainingSets(pathSets.begin() + 1, pathSets.end());

            pathSets.clear();

            for (const auto &set : remainingSets)
            {
                // Create the disjoint sets from the selectedSet and each set in the remaining sets
                DisjointSets disjointSets = makeDisjointSet(selectedSet, set);

                // Insert the disjoint sets to the pathSets
                pathSets.insert(pathSets.end(), disjointSets.begin(), disjointSets.end());
            }
        }

        return probaSets;
    }

    DebugInfo toProbaSetDebug(NodeID src, NodeID dst, PathSets pathSets)
    {   
        // Initialize the debug information
        DebugInfo debugInfo;

        // Initialize the probability sets
        ProbaSets probaSets;

        // Reserve space for the probability sets
        probaSets.reserve(pathSets.size() * 3); // Assuming each path set can generate up to 3 disjoint sets

        // Initialize the iteration counter
        int iteration = 0;

        // convert the path sets to the probability sets
        while (pathSets.size() > 0)
        {   

            // Measure the time taken for each iteration
            auto start = std::chrono::high_resolution_clock::now();

            if (pathSets.size() == 1)
            {
                probaSets.push_back(pathSets[0]);
                break;
            }

            Set selectedSet = pathSets.front();

            probaSets.push_back(selectedSet);

            std::vector<std::vector<int>> remainingSets(pathSets.begin() + 1, pathSets.end());

            pathSets.clear();

            for (const auto &set : remainingSets)
            {
                // Create the disjoint sets from the selectedSet and each set in the remaining sets
                DisjointSets disjointSets = makeDisjointSet(selectedSet, set);

                // Insert the disjoint sets to the pathSets
                pathSets.insert(pathSets.end(), disjointSets.begin(), disjointSets.end());
            }

            // Measure the end time of the iteration
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

    double probaSetToAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const ProbaSets &probaSets)
    {
        // Save the final result
        double avail = 0.0;
        for (const auto &set : probaSets)
        {
            // Save the temp result
            double temp = 1.0;
            for (const auto &num : set)
            {
                temp *= probaMap[num];
            }
            avail += temp;
        }

        // return the availability
        return avail;
    }

    double evalAvail(NodeID src, NodeID dst, const ProbabilityMap &probaMap, const PathSets &pathSets)
    {
        // Convert the path sets to the probability sets
        ProbaSets probaSets = toProbaSet(src, dst, pathSets);

        // Compute the availability
        return probaSetToAvail(src, dst, probaMap, probaSets);
    }

    std::vector<AvailTriple> evalAvailTopo(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<PathSets> &pathsetsList)
    {
        std::vector<AvailTriple> availList;

        for (size_t i = 0; i < nodePairs.size(); ++i)
        {   
            const auto &[src, dst] = nodePairs[i];
            const auto &pathsets = pathsetsList[i];

            double availability = evalAvail(src, dst, probaMap, pathsets);
            availList.emplace_back(src, dst, availability);
        }

        return availList;
    }

    std::vector<AvailTriple> evalAvailTopoParallel(const NodePairs &nodePairs, const ProbabilityMap &probaMap, const std::vector<PathSets> &pathsetsList)
    {
        std::vector<AvailTriple> availList(nodePairs.size());

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < nodePairs.size(); ++i)
        {
            const auto &[src, dst] = nodePairs[i];
            const auto &pathsets = pathsetsList[i];

            double availability = evalAvail(src, dst, probaMap, pathsets);
            availList[i] = std::make_tuple(src, dst, availability);
        }

        return availList;
    }

} // namespace pyrbdpp::pathset