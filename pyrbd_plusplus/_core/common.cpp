#include <pyrbd_plusplus/common.hpp>
#include <algorithm>

namespace pyrbdpp
{
    DisjointSets makeDisjointSet(const Set &set1, Set set2)
    {
        // Check disjointness of set1 and set2: if x in set1 and -x in set2, set1 and set2 are disjoint
        // RC set: the elements in set1 but not in set2
        Set RC;

        // Check if set1 and set2 are disjoint
        for (const int &elem : set1)
        {
            if (std::find(set2.begin(), set2.end(), -elem) != set2.end())
            {
                return {set2};
            }
            if (std::find(set2.begin(), set2.end(), elem) == set2.end())
            {
                RC.push_back(elem);
            }
        }

        // Check if RC is empty
        if (RC.empty())
        {
            return DisjointSets{};
        }

        // Save the result
        DisjointSets result;

        // Add the elements in RC to set2 to create several new sets
        for (size_t i = 0; i < RC.size(); i++)
        {
            // every time add one more element from RC to set2 and the element to add is negative.
            set2.push_back(-RC[i]);
            // save the new set
            result.push_back(set2);
            // inverse the last element to keep only the last element negative
            set2.back() = -set2.back();
        }

        return result;
    }
} // namespace pyrbdpp