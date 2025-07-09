#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

namespace pyrbdpp
{   
    // Declaration of the short types
    using NodeID = int;
    using Set = std::vector<NodeID>;
    using DisjointSets = std::vector<Set>;
    using ProbaSets    = std::vector<Set>;
    using DebugInfo   = std::map<NodeID, std::pair<NodeID,double>>;
    using AvailTriple = std::tuple<NodeID,NodeID,double>;
    using NodePairs   = std::vector<std::pair<NodeID,NodeID>>;


    // Declaration of the Probability Map class.
    class ProbabilityMap
    {

    private:
        // Two arrays to store the probabilities for positive and negative indices.
        std::vector<double> pos_array;
        std::vector<double> neg_array;

    public:
        ProbabilityMap(const std::map<int, double> avail_arr) : pos_array(avail_arr.size() + 1), neg_array(avail_arr.size() + 1)
        {
            for (const auto &pair : avail_arr)
            {
                pos_array[pair.first] = pair.second;
                neg_array[pair.first] = 1 - pair.second;
            }
        }

        ProbabilityMap(std::initializer_list<std::pair<const int, double>> init_list)
        {
            int max_id = 0;
            for (const auto& pair : init_list) {
                max_id = std::max(max_id, pair.first);
            }
            
            pos_array.resize(max_id + 1);
            neg_array.resize(max_id + 1);

            for (const auto &pair : init_list)
            {
                pos_array[pair.first] = pair.second;
                neg_array[pair.first] = 1 - pair.second;
            }
        }

        double operator[](int i) const
        {
            if (std::abs(i) > pos_array.size())
            {
                throw std::out_of_range("Index out of range in ProbabilityArray");
            }
            if (i == 0)
            {
                std::cerr << "Error on handling node 0!" << std::endl;
            }
            if (i > 0)
            {
                return pos_array[i];
            }
            else
            {
                return neg_array[-i];
            }
        }

        void print() const
        {
            std::cout << "Positive Array: ";
            for (double elem : pos_array)
            {

                std::cout << elem << " ";
            }
            std::cout << std::endl;

            std::cout << "Negative Array: ";
            for (double elem : neg_array)
            {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
    };


    // Declaration of the SDP class
    class SDP
    {
    private:
        bool is_complementary;
        std::vector<int> set;

    public:
        SDP(bool is_comp, std::vector<int> s)
            : is_complementary(is_comp), set(std::move(s)) {}

        SDP() : is_complementary(false), set() {}

        bool isComplementary() const { return is_complementary; }

        const std::vector<int> &getSet() const { return set; }

        std::vector<int>::iterator begin() { return set.begin(); }

        std::vector<int>::iterator end() { return set.end(); }

        std::vector<int>::const_iterator begin() const { return set.begin(); }

        std::vector<int>::const_iterator end() const { return set.end(); }

        int operator[](size_t idx) const { return set[idx]; }

        size_t size() const { return set.size(); }

        void remove(int elem)
        {
            set.erase(std::remove(set.begin(), set.end(), elem), set.end());
        }

        bool equals(const SDP &other) const
        {
            if (is_complementary != other.is_complementary || set.size() != other.set.size())
            {
                return false;
            }
            return std::equal(set.begin(), set.end(), other.set.begin());
        }

        friend std::ostream &operator<<(std::ostream &os, const SDP &sdp)
        {
            if (sdp.is_complementary)
            {
                os << "-";
            }
            os << "{ ";
            for (const auto &elem : sdp.set)
            {
                os << elem << " ";
            }
            os << "}";
            return os;
        }

        // Print the SDP set
        void print() const
        {
            if (is_complementary)
            {
                std::cout << "-";
            }
            std::cout << "{ ";
            for (const auto &elem : set)
            {
                std::cout << elem << " ";
            }
            std::cout << "}";
        }
    };

    /**
     * @brief Create a couple of disjoint sets from set2
     * Algorithm:
     * 1. Find the set: set1 \ set2, e.g. set1 = {1, 2, 3, 4, 5}, set2 = {2, 4}, then RC = set1 \ set2 = {1, 3, 5}
     * 2. add the elements in RC to set2 to create several new sets. The rule is every time add one more element from RC to set2 and the last element to add is negative.
     *      1) first step: {2, 4} -> {2, 4, -1}
     *      2) second step: {2, 4, 1} -> {2, 4, 1, -3}
     *      3) third step: {2, 4, 1, -3} -> {2, 4, 1, -3, -5}
     * 3. return the new disjoint sets {{2, 4, 1}, {2, 4, 1, -3}, {2, 4, 1, 3, -5}}
     * @param set1 First set
     * @param set2 Second set
     * @return Disjointed sets
     *
     */
    DisjointSets makeDisjointSet(const Set &set1, Set set2);

}