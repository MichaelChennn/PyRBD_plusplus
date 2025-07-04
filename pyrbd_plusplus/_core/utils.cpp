#include <pyrbd_plusplus/utils.hpp>
#include <algorithm>

namespace pyrbdpp::utils
{
   
    bool isSubSet(const Set &set1, const Set &set2)
    {
        // Check if set1 is a subset of set2
        for (const auto &elem : set1)
        {
            if (std::find(set2.begin(), set2.end(), elem) == set2.end())
            {
                return false;
            }
        }
        return true;
    }
}